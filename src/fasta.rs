use crossbeam::atomic::AtomicCell;

use std::sync::{Arc, RwLock};

use std::thread;
use std::thread::Builder;
use std::thread::JoinHandle;

use std::fs;
use std::fs::File;
use std::io::{BufReader, Read, BufRead};

use crossbeam::queue::{ArrayQueue, PushError};
use crossbeam::utils::Backoff;

const STACKSIZE: usize = 256 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)

#[derive(PartialEq)]
pub struct Sequence {
    pub seq: Vec<u8>,
    pub id:  String,
}

#[derive(PartialEq)]
pub enum ThreadCommand<T> {
    Work(T),
    Terminate,
}

impl ThreadCommand<Sequence> {
    // Consumes the ThreadCommand, which is just fine...
    pub fn unwrap(self) -> Sequence {
        match self {
            ThreadCommand::Work(x)   => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
}

pub fn shuffle_fasta_file(filename: &str, threads: usize, buckets: usize) {

    // Start this first, just so everything is ready to go!
    let (seq_queue, seqs_submitted, generator_done, generator_joinhandle) = sequence_generator(filename, threads*32);

    match fs::create_dir("shuffle_round1") {
        Ok(_)   => (),
        Err(y) => panic!("Following error when trying to create shuffle_round1 directory: {}", y)
    };

    match fs::create_dir("shuffle_round2") {
        Ok(_)   => (),
        Err(y) => panic!("Following error when trying to create shuffle_round2 directory: {}", y)
    };

    match fs::create_dir("shuffle_round3") {
        Ok(_)   => (),
        Err(y) => panic!("Following error when trying to create shuffle_round3 directory: {}", y)
    };

    match fs::create_dir("shuffle_round4") {
        Ok(_)   => (),
        Err(y) => panic!("Following error when trying to create shuffle_round4 directory: {}", y)
    };

    match fs::create_dir("shuffle_round5") {
        Ok(_)   => (),
        Err(y) => panic!("Following error when trying to create shuffle_round5 directory: {}", y)
    };



}


// Takes a file and submits Sequence type to a buffer...
// Fills up the seq_buffer that it returns...
// Will park if output buffer is full
pub fn sequence_generator(
    filename: &str,
    queue_size: usize,
) -> (Arc<ArrayQueue<ThreadCommand<Sequence>>>, 
      Arc<AtomicCell<usize>>, 
      Arc<RwLock<bool>>, 
      JoinHandle<()>, )
{
    let seqs_processed = Arc::new(AtomicCell::new(0 as usize));
    let seq_queue = Arc::new(ArrayQueue::<ThreadCommand<Sequence>>::new(queue_size));
    let generator_done = Arc::new(RwLock::new(false));

    let generator;

    let filename = filename.to_string();

    { // Explicit lifetime
        let generator_done = Arc::clone(&generator_done);
        let seq_queue = Arc::clone(&seq_queue);
        let seqs_processed = Arc::clone(&seqs_processed);

        generator = thread::Builder::new()
                            .name("Generator".to_string())
                            .stack_size(STACKSIZE)
                            .spawn(move||
        {
            let mut id: String = String::from("INVALID_ID_FIRST_ENTRY_YOU_SHOULD_NOT_SEE_THIS");
            let mut seqbuffer: Vec<u8> = Vec::with_capacity(8 * 1024 * 1024); // 8 Mb to start, will likely increase...
            let mut seqlen: usize = 0;
            let mut buffer: Vec<u8> = Vec::with_capacity(8192*10);

            let file = match File::open(&filename) {
                Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
                Ok(file) => file,
            };

            let file = BufReader::with_capacity(32 * 1024 * 1024, file);

            let fasta: Box<dyn Read> = 
                if filename.ends_with("gz") {
                    Box::new(flate2::read::GzDecoder::new(file))
                } else if filename.ends_with("snappy") {
                    Box::new(snap::Reader::new(file))
                } else {
                    Box::new(file)
                };

            let mut reader = BufReader::with_capacity(32 * 1024 * 1024, fasta);

            let backoff = Backoff::new();

            while let Ok(bytes_read) = reader.read_until(b'\n', &mut buffer) {

                if bytes_read == 0 { // No more reads, thus no more data...
                    // Submit the last sequence (or in the case of some genomes, the entire sequence)
                    seqs_processed.fetch_add(1 as usize);
                    let wp = ThreadCommand::Work(Sequence { seq: seqbuffer[..seqlen].to_vec(), id: id });
                    seqbuffer.clear();

                    let mut result = seq_queue.push(wp);
                    while let Err(PushError(wp)) = result {
                        if backoff.is_completed() {
                            thread::park();
                        } else {
                            backoff.snooze();
                        }
                        result = seq_queue.push(wp);
                    }
                    break;
                }

                match buffer[0] {
                    // 62 is a > meaning we have a new sequence id.
                    62 => {
                        seqs_processed.fetch_add(1 as usize);
                        let wp = ThreadCommand::Work(Sequence { seq: seqbuffer[..seqlen].to_vec(), id: id });
                        seqbuffer.clear();
                        seqlen = 0;

                        let mut result = seq_queue.push(wp);
                        while let Err(PushError(wp)) = result {
                            // TODO: Same code, two spots... Make into a proper function
                            if backoff.is_completed() {
                                thread::park();
                            } else {
                                backoff.snooze();
                            }
                            result = seq_queue.push(wp);
                        }

                        let slice_end = bytes_read.saturating_sub(1);
                        id = String::from_utf8(buffer[1..slice_end].to_vec()).expect("Invalid UTF-8 encoding...");
                    },
                    _  => {
                        let slice_end = bytes_read.saturating_sub(1);
                        seqbuffer.extend_from_slice(&buffer[0..slice_end]);
                        seqlen = seqlen.saturating_add(slice_end);
                    }
                }
            }
            // buffer.clear();
            *generator_done.write().unwrap() = true;
        }).unwrap();
    }
    (seq_queue, seqs_processed, generator_done, generator)
}
