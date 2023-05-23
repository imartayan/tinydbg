use seq_io::fasta::{Reader, RefRecord};
use seq_io::parallel::read_process_fasta_records;
use seq_io::BaseRecord;
use std::fs::File;
use std::path::Path;
use std::slice::Iter;

pub trait ReadProcess {
    fn process<F: FnMut(Iter<u8>)>(self, f: F);
    fn parallel_process<F: Send + Sync + Fn(Iter<u8>)>(self, threads: u32, queue_len: usize, f: F);
}

pub struct Fasta {
    reader: Reader<File>,
}

impl Fasta {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Self {
        Self {
            reader: Reader::from_path(path).expect("Failed to open file"),
        }
    }
}

impl ReadProcess for Fasta {
    fn process<F: FnMut(Iter<u8>)>(mut self, mut f: F) {
        while let Some(result) = self.reader.next() {
            let record = result.expect("Error reading record");
            f(record.seq().iter());
        }
    }

    fn parallel_process<F: Send + Sync + Fn(Iter<u8>)>(self, threads: u32, queue_len: usize, f: F) {
        read_process_fasta_records(
            self.reader,
            threads,
            queue_len,
            |record: RefRecord, _: &mut Option<()>| {
                f(record.seq().iter());
            },
            |_, _| None::<()>,
        )
        .unwrap();
    }
}
