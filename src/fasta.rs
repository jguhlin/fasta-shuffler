use crossbeam::atomic::AtomicCell;

use std::sync::{Arc, RwLock};

use std::thread;
use std::thread::Builder;
use std::thread::JoinHandle;

use std::fs::File;
use std::io::{BufReader, Read, BufRead};

use crossbeam::queue::{ArrayQueue, PushError};
use crossbeam::utils::Backoff;

#[derive(PartialEq)]
pub struct Sequence {
    pub rawseq: Vec<u8>,
    pub id:     String,
}

#[derive(PartialEq)]
pub enum ThreadCommand<T> {
    Work(T),
    Terminate,
}

#[derive(PartialEq)]
pub struct SequenceTargetContexts {
    pub id: String,
    pub target: String,
    pub contexts: Vec<String>
}

#[derive(PartialEq)]
pub struct SequenceKmers {
    pub id: String,
    pub kmers: Vec<String>
}

// TODO: Should be able to combine these with an ENUM
pub type SequenceBatch = Vec<SequenceTargetContexts>;
pub type SequenceBatchKmers = Vec<SequenceKmers>;

impl ThreadCommand<Sequence> {
    // Consumes the ThreadCommand, which is just fine...
    pub fn unwrap(self) -> Sequence {
        match self {
            ThreadCommand::Work(x)   => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
}

impl ThreadCommand<SequenceTargetContexts> {
    // Consumes the ThreadCommand, which is just fine...
    pub fn unwrap(self) -> SequenceTargetContexts {
        match self {
            ThreadCommand::Work(x)   => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
}

impl ThreadCommand<SequenceKmers> {
    // Consumes the ThreadCommand, which is just fine...
    pub fn unwrap(self) -> SequenceKmers {
        match self {
            ThreadCommand::Work(x)   => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
}

impl ThreadCommand<SequenceBatch> {
    // Consumes the ThreadCommand, which is just fine...
    pub fn unwrap(self) -> SequenceBatch {
        match self {
            ThreadCommand::Work(x)   => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
}

impl ThreadCommand<SequenceBatchKmers> {
    // Consumes the ThreadCommand, which is just fine...
    pub fn unwrap(self) -> SequenceBatchKmers {
        match self {
            ThreadCommand::Work(x)   => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
}

// Takes a file and submits Sequence type to buffers...
// Fills up the seq_buffer that it returns
// UP TO the calling fn/thread to handle clean-up...
pub fn sequence_generator(
    kmer_size: usize,
    filename: &str,
    epochs: u64
) -> (Arc<ArrayQueue<ThreadCommand<Sequence>>>, Arc<AtomicCell<usize>>, Arc<RwLock<bool>>, JoinHandle<()>, Vec<JoinHandle<()>>)
// seq_queue jobs generator_done generator children

{
    let jobs = Arc::new(AtomicCell::new(0 as usize));
    let seq_queue = Arc::new(ArrayQueue::<ThreadCommand<Sequence>>::new(256));
    let rawseq_queue = Arc::new(ArrayQueue::<ThreadCommand<Sequence>>::new(256));
    let generator_done = Arc::new(RwLock::new(false));

    let generator;

    let mut children = Vec::new();

    let filename = filename.to_string();

    // IO-bound, so 4 threads never seems to hurt anything (but seems to help)
    for _ in 0..4 {
        let seq_queue = Arc::clone(&seq_queue);
        let rawseq_queue = Arc::clone(&rawseq_queue);
        let jobs = Arc::clone(&jobs);

        let child = match Builder::new()
                        .name("IOWorker".into())
                        .stack_size(WORKERSTACKSIZE)
                        .spawn(move || io_worker_thread(kmer_size, rawseq_queue, seq_queue, jobs)) {
                            Ok(x)  => x,
                            Err(y) => panic!("{}", y)
                        };
        
        children.push(child);
    }

    { // Explicit lifetime
        let generator_done = Arc::clone(&generator_done);
        let rawseq_queue = Arc::clone(&rawseq_queue);
        let jobs = Arc::clone(&jobs);
        let mut buffer: Vec<u8> = Vec::with_capacity(1024);

        generator = thread::Builder::new()
                            .name("Generator".to_string())
                            .stack_size(STACKSIZE)
                            .spawn(move||
        {
            for _ in 0..epochs {
                let mut id: String = String::from("INVALID_ID_FIRST_ENTRY_YOU_SHOULD_NOT_SEE_THIS");
                let mut seqbuffer: Vec<u8> = Vec::with_capacity(8 * 1024 * 1024); // 8 Mb to start, will likely increase...
                let mut seqlen: usize = 0;

                let file = match File::open(&filename) {
                    Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
                    Ok(file) => file,
                };

                let file = BufReader::with_capacity(64 * 1024 * 1024, file);

                let fasta: Box<dyn Read> = 
                    if filename.ends_with("gz") {
                        Box::new(flate2::read::GzDecoder::new(file))
            		} else if filename.ends_with("snappy") {
		                Box::new(snap::Reader::new(file))
                    } else {
                        Box::new(file)
                    };

                let mut reader = BufReader::with_capacity(128 * 1024 * 1024, fasta);

                let backoff = Backoff::new();

                while let Ok(bytes_read) = reader.read_until(b'\n', &mut buffer) {

                    if bytes_read == 0 { // No more reads, thus no more data...
                        // Submit the last sequence (or in the case of some genomes, the entire sequence)
                        jobs.fetch_add(1 as usize);
                        let wp = ThreadCommand::Work(Sequence { rawseq: seqbuffer[..seqlen].to_vec(), id: id });
                        seqbuffer.clear();

                        let mut result = rawseq_queue.push(wp);
                        while let Err(PushError(wp)) = result {
                            result = rawseq_queue.push(wp);
                        }

                        backoff.spin();
                        backoff.spin();
                        backoff.spin();
                        backoff.spin(); // Very slight delay then issue terminate commands...

                        break;
                    }

                    match buffer[0] {
                        // 62 is a > meaning we have a new sequence id.
                        62 => {
                            jobs.fetch_add(1 as usize);
                            let wp = ThreadCommand::Work(Sequence { rawseq: seqbuffer[..seqlen].to_vec(), id: id });
                            seqbuffer.clear();
                            seqlen = 0;

                            let mut result = rawseq_queue.push(wp);
                            while let Err(PushError(wp)) = result {
                                result = rawseq_queue.push(wp);
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

                buffer.clear();
                }

            }
            *generator_done.write().unwrap() = true;
            for _ in 0..4 {
                let mut result = rawseq_queue.push(ThreadCommand::Terminate);
                while let Err(PushError(wp)) = result {
                    // println!("Parking generator thread...");
                    thread::park();
                    result = rawseq_queue.push(wp);
                }
            }}
            
            ).unwrap();
    }


    (seq_queue, jobs, generator_done, generator, children)
}

fn io_worker_thread(
    kmer_size: usize,
    rawseq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    jobs: Arc<AtomicCell<usize>>) {

    let backoff = Backoff::new();

    loop {
        if let Ok(command) = rawseq_queue.pop() {
            // We got the work packet!

            // We are finished, end the thread...
            if let ThreadCommand::Terminate = command {
                return;
            }

            let seqpacket = command.unwrap();
            jobs.fetch_sub(1);
            let mut rawseq = seqpacket.rawseq;
            let id         = seqpacket.id;

            
            for (start_coords, end_coords) in coords {

                if (end_coords - start_coords) >= kmer_size {

                    jobs.fetch_add(1 as usize);
                    let wp = ThreadCommand::Work(Sequence { rawseq: rawseq[start_coords..end_coords].to_vec(), id: id.clone() });

                    let mut result = seq_queue.push(wp);
                    while let Err(PushError(wp)) = result {
                        backoff.snooze();
                        result = seq_queue.push(wp);
                    }

                    let mut rc = rawseq[start_coords..end_coords].to_vec();
                    complement_nucleotides(&mut rc);
                    rc.reverse();

                    jobs.fetch_add(1 as usize);
                    let wp = ThreadCommand::Work(Sequence { rawseq: rc.clone(), id: id.clone() });

                    let mut result = seq_queue.push(wp);
                    while let Err(PushError(wp)) = result {
                        // println!("Parking IO worker thread -- FULL");
                        thread::park();
                        result = seq_queue.push(wp);
                    }


                }
                
            }

            backoff.reset();

        } else {
            backoff.snooze();
        }
    }
}
