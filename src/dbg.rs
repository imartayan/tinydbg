use crate::kmer::{Base, IntKmer, Kmer};
use ahash::{HashSet, HashSetExt};
use anyhow::Result;
use std::collections::BTreeSet;
use std::io::{Read, Write};
use std::marker::PhantomData;
use sucds::bit_vectors::BitVector;
use sucds::mii_sequences::{EliasFano, EliasFanoBuilder};
use sucds::Serializable;

pub trait Dbg<const K: usize, T, KT>
where
    T: Base,
    KT: Kmer<K, T>,
{
    fn contains(&self, kmer: KT) -> bool;
    fn successors(&self, kmer: KT) -> Vec<KT> {
        kmer.successors()
            .into_iter()
            .filter(|&s| self.contains(s))
            .collect()
    }
}

pub trait DbgBuilder<const K: usize, T, KT, DT>
where
    T: Base,
    KT: Kmer<K, T>,
    DT: Dbg<K, T, KT>,
{
    fn new() -> Self;
    fn insert(self, kmer: IntKmer<K, T>) -> Self;
    fn build(self) -> DT;
}

pub struct HashDbg<const K: usize, T: Base> {
    data: HashSet<T>,
}

pub type HashDbgBuilder<const K: usize, T> = HashDbg<K, T>;

pub struct DenseDbg<const K: usize, T: Base> {
    data: BitVector,
    phantom: PhantomData<T>,
}

pub type DenseDbgBuilder<const K: usize, T> = DenseDbg<K, T>;

pub struct SparseDbg<const K: usize, T: Base> {
    data: EliasFano,
    phantom: PhantomData<T>,
}

pub struct SparseDbgBuilder<const K: usize, T: Base> {
    positions: BTreeSet<T>,
}

macro_rules! impl_traits {
($($t:ty),+) => {$(
    impl<const K: usize> Dbg<K, $t, IntKmer<K, $t>> for HashDbg<K, $t> {
        fn contains(&self, kmer: IntKmer<K, $t>) -> bool {
            self.data.contains(&kmer.to_int())
        }
    }

    impl<const K: usize> DbgBuilder<K, $t, IntKmer<K, $t>, HashDbg<K, $t>> for HashDbgBuilder<K, $t> {
        fn new() -> Self {
            Self {
                data: HashSet::new(),
            }
        }

        fn insert(mut self, kmer: IntKmer<K, $t>) -> Self {
            self.data.insert(kmer.to_int());
            self
        }

        fn build(self) -> HashDbg<K, $t> {
            self
        }
    }

    impl<const K: usize> Dbg<K, $t, IntKmer<K, $t>> for DenseDbg<K, $t> {
        fn contains(&self, kmer: IntKmer<K, $t>) -> bool {
            let pos = kmer.to_int() as usize;
            self.data.get_bit(pos).expect("Out of bounds")
        }
    }

    impl<const K: usize> DbgBuilder<K, $t, IntKmer<K, $t>, DenseDbg<K, $t>> for DenseDbgBuilder<K, $t> {
        fn new() -> Self {
            Self {
                data: BitVector::from_bit(false, 1 << (2 * K)),
                phantom: PhantomData,
            }
        }

        fn insert(mut self, kmer: IntKmer<K, $t>) -> Self {
            let pos = kmer.to_int() as usize;
            self.data.set_bit(pos, true).expect("Out of bounds");
            self
        }

        fn build(self) -> DenseDbg<K, $t> {
            self
        }
    }

    impl<const K: usize> Dbg<K, $t, IntKmer<K, $t>> for SparseDbg<K, $t> {
        fn contains(&self, kmer: IntKmer<K, $t>) -> bool {
            let pos = kmer.to_int() as usize;
            if let Some(rank) = self.data.rank(pos) {
                self.data.select(rank) == Some(pos)
            } else {
                false
            }
        }
    }

    impl<const K: usize> DbgBuilder<K, $t, IntKmer<K, $t>, SparseDbg<K, $t>> for SparseDbgBuilder<K, $t> {
        fn new() -> Self {
            Self {
                positions: BTreeSet::new(),
            }
        }

        fn insert(mut self, kmer: IntKmer<K, $t>) -> Self {
            self.positions.insert(kmer.to_int());
            self
        }

        fn build(self) -> SparseDbg<K, $t> {
            let mut efb = EliasFanoBuilder::new(1 << (2 * K), self.positions.len())
                .expect("Failed to create Elias-Fano Builder");
            efb.extend(self.positions.iter().map(|&x| x as usize)).unwrap();
            SparseDbg {
                data: efb.build().enable_rank(),
                phantom: PhantomData,
            }
        }
    }
)*}}

impl_traits!(u8, u16, u32, u64, u128);

impl<const K: usize, T: Base> Serializable for DenseDbg<K, T> {
    fn serialize_into<W: Write>(&self, mut writer: W) -> Result<usize> {
        self.data.serialize_into(&mut writer)
    }

    fn deserialize_from<R: Read>(mut reader: R) -> Result<Self> {
        Ok(Self {
            data: BitVector::deserialize_from(&mut reader)?,
            phantom: PhantomData,
        })
    }

    fn size_in_bytes(&self) -> usize {
        self.data.size_in_bytes()
    }
}

impl<const K: usize, T: Base> Serializable for SparseDbg<K, T> {
    fn serialize_into<W: Write>(&self, mut writer: W) -> Result<usize> {
        self.data.serialize_into(&mut writer)
    }

    fn deserialize_from<R: Read>(mut reader: R) -> Result<Self> {
        Ok(Self {
            data: EliasFano::deserialize_from(&mut reader)?,
            phantom: PhantomData,
        })
    }

    fn size_in_bytes(&self) -> usize {
        self.data.size_in_bytes()
    }
}
