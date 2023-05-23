use std::env;
use tinydbg::dbg::{Dbg, DbgBuilder, HashDbgBuilder};
use tinydbg::kmer::{Canonical, IntKmer, Kmer};
use tinydbg::reads::{Fasta, ReadProcess};

const K: usize = 31;
type T = u64;

fn main() {
    let args: Vec<String> = env::args().collect();
    let filename = args.get(1).expect("No filename given").as_str();
    let reads = Fasta::from_file(filename);
    let mut builder = HashDbgBuilder::<K, T>::new();
    reads.process(|bytes| {
        IntKmer::<K, T>::iter_from_chars(bytes).for_each(|kmer| builder.insert(kmer.canonical()));
    });
    let dbg = builder.build();
    let reads = Fasta::from_file(filename);
    reads.process(|bytes| {
        IntKmer::<K, T>::iter_from_chars(bytes)
            .for_each(|kmer| assert!(dbg.contains(kmer.canonical())));
    });
}
