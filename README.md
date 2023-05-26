# tinydbg

## What is tinydbg?

It's a Rust library providing tools to build succinct De Bruijn graphs.
It uses [sucds](https://docs.rs/sucds/) for the underlying succinct data structures (bit vector and Elias-Fano encoding).

Currently, it provides 3 data structures for De Bruijn graphs:
- a "dense" structure that uses a bit vector; it is mostly suited for small k-mers (k <= 16) since it allocates an array of 4^k bits
- a "sparse" structure that uses Elias-Fano encoding; it is mostly suited for long k-mers, when most k-mers are not expected to be contained in the graph
- a classic structure using a HashSet (from [ahash](https://docs.rs/ahash/)); this one is not particularly optimized and is mostly used for comparison with the other structures

## Structure

There are currently 3 modules:
- `kmer` provides an interface for manipulating k-mers; for better performance, k-mers are represented as raw integers under the hood
- `dbg` implements the data structures for De Bruijn graphs presented above
- `reads` is a wrapper around [seq_io](https://docs.rs/seq_io/0.4.0-alpha.0/) for parsing fasta files and processing them easily

## Examples

Small examples are available in the `examples` folder.
They illustrate how to insert every k-mer of a fasta file in dense / sparse / hash-based structure, and how to verify that they are actually contained in the structure.

For instance, the dense structure can be tested with:
```sh
cargo r -r --example dense -- <sequence.fasta>
```

## Reference

[Conway & Bromage - Succinct data structures for assembling large genomes](https://academic.oup.com/bioinformatics/article/27/4/479/198367)
