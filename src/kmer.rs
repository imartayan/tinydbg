use num_traits::int::PrimInt;
use std::{iter::FilterMap, marker::PhantomData};

pub trait Base: PrimInt {
    fn from_char(b: u8) -> Option<Self>;
    fn to_char(self) -> u8;
    fn bases() -> [Self; 4];
}

pub trait Kmer<const K: usize, T: Base>: Sized + Copy {
    fn from_int(s: T) -> Self;
    fn to_int(self) -> T;
    fn empty() -> Self;
    fn mask() -> T;
    fn extend(self, base: T) -> Self;
    fn append(self, base: T) -> Self;
    #[inline]
    fn successors(self) -> [Self; 4] {
        T::bases().map(|base| self.append(base))
    }
    fn from_chars(bytes: &[u8]) -> Self {
        bytes
            .iter()
            .filter_map(|&b| T::from_char(b))
            .take(K)
            .fold(Self::empty(), |s, base| s.extend(base))
    }
    fn to_chars(self) -> [u8; K];
    fn iter_from_bases<I: Iterator<Item = T>>(bases: I) -> KmerIterator<K, T, I, Self> {
        KmerIterator {
            kmer: Self::empty(),
            bases,
            init: false,
            phantom: PhantomData,
        }
    }
    fn iter_from_chars<I: Iterator<Item = u8>>(
        bytes: I,
    ) -> KmerIterator<K, T, FilterMap<I, fn(u8) -> Option<T>>, Self> {
        Self::iter_from_bases(bytes.filter_map(T::from_char))
    }
}

pub trait Canonical<const K: usize, T: Base>: Kmer<K, T> {
    fn rev_comp(self) -> Self;
    fn canonical(self) -> Self {
        let rc = self.rev_comp();
        if self.to_int() < rc.to_int() {
            self
        } else {
            rc
        }
    }
}

pub struct KmerIterator<const K: usize, T, I, KT>
where
    T: Base,
    I: Iterator<Item = T>,
    KT: Kmer<K, T>,
{
    kmer: KT,
    bases: I,
    init: bool,
    phantom: PhantomData<T>,
}

impl<const K: usize, T, I, KT> Iterator for KmerIterator<K, T, I, KT>
where
    T: Base,
    I: Iterator<Item = T>,
    KT: Kmer<K, T>,
{
    type Item = KT;
    fn next(&mut self) -> Option<Self::Item> {
        if !self.init {
            self.init = true;
            for _ in 0..K {
                if let Some(base) = self.bases.next() {
                    self.kmer = self.kmer.extend(base);
                } else {
                    return None;
                }
            }
            Some(self.kmer)
        } else {
            if let Some(base) = self.bases.next() {
                self.kmer = self.kmer.append(base);
                Some(self.kmer)
            } else {
                None
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct IntKmer<const K: usize, T: Base>(T);

macro_rules! impl_traits {
($($t:ty),+) => {$(
    impl Base for $t {
        #[inline]
        fn from_char(b: u8) -> Option<Self> {
            match b {
                b'A' => Some(0b00),
                b'C' => Some(0b01),
                b'G' => Some(0b10),
                b'T' => Some(0b11),
                _ => None,
            }
        }
        #[inline]
        fn to_char(self) -> u8 {
            match self {
                0b00 => b'A',
                0b01 => b'C',
                0b10 => b'G',
                _ => b'T',
            }
        }
        #[inline]
        fn bases() -> [Self; 4] {
            [0, 1, 2, 3]
        }
    }

    impl<const K: usize> Kmer<K, $t> for IntKmer<K, $t> {
        #[inline(always)]
        fn from_int(s: $t) -> Self {
            Self(s)
        }
        #[inline(always)]
        fn to_int(self) -> $t {
            self.0
        }
        #[inline(always)]
        fn empty() -> Self {
            Self::from_int(0)
        }
        #[inline(always)]
        fn mask() -> $t {
            (1 << (2 * K)) - 1
        }
        #[inline]
        fn extend(self, base: $t) -> Self {
            Self::from_int((self.to_int() << 2) | base)
        }
        #[inline]
        fn append(self, base: $t) -> Self {
            Self::from_int(((self.to_int() << 2) | base) & Self::mask())
        }
        fn to_chars(self) -> [u8; K] {
            let mut res = [0u8; K];
            let mut s = self.to_int();
            for i in 0..K {
                res[K - i - 1] = (s & 0b11).to_char();
                s >>= 2;
            }
            res
        }
    }
)*}}

impl_traits!(u8, u16, u32, u64, u128);

impl<const K: usize> Canonical<K, u8> for IntKmer<K, u8> {
    fn rev_comp(self) -> Self {
        let mut res = !self.to_int();
        res = (res >> 2 & 0x33) | (res & 0x33) << 2;
        res = (res >> 4 & 0x0F) | (res & 0x0F) << 4;
        Self::from_int(res >> (2 * (4 - K)))
    }
}

impl<const K: usize> Canonical<K, u16> for IntKmer<K, u16> {
    fn rev_comp(self) -> Self {
        let mut res = !self.to_int();
        res = (res >> 2 & 0x3333) | (res & 0x3333) << 2;
        res = (res >> 4 & 0x0F0F) | (res & 0x0F0F) << 4;
        res = (res >> 8 & 0x00FF) | (res & 0x00FF) << 8;
        Self::from_int(res >> (2 * (8 - K)))
    }
}

impl<const K: usize> Canonical<K, u32> for IntKmer<K, u32> {
    fn rev_comp(self) -> Self {
        let mut res = !self.to_int();
        res = (res >> 2 & 0x3333_3333) | (res & 0x3333_3333) << 2;
        res = (res >> 4 & 0x0F0F_0F0F) | (res & 0x0F0F_0F0F) << 4;
        res = (res >> 8 & 0x00FF_00FF) | (res & 0x00FF_00FF) << 8;
        res = (res >> 16 & 0x0000_FFFF) | (res & 0x0000_FFFF) << 16;
        Self::from_int(res >> (2 * (16 - K)))
    }
}

impl<const K: usize> Canonical<K, u64> for IntKmer<K, u64> {
    fn rev_comp(self) -> Self {
        let mut res = !self.to_int();
        res = (res >> 2 & 0x3333_3333_3333_3333) | (res & 0x3333_3333_3333_3333) << 2;
        res = (res >> 4 & 0x0F0F_0F0F_0F0F_0F0F) | (res & 0x0F0F_0F0F_0F0F_0F0F) << 4;
        res = (res >> 8 & 0x00FF_00FF_00FF_00FF) | (res & 0x00FF_00FF_00FF_00FF) << 8;
        res = (res >> 16 & 0x0000_FFFF_0000_FFFF) | (res & 0x0000_FFFF_0000_FFFF) << 16;
        res = (res >> 32 & 0x0000_0000_FFFF_FFFF) | (res & 0x0000_0000_FFFF_FFFF) << 32;
        Self::from_int(res >> (2 * (32 - K)))
    }
}

impl<const K: usize> Canonical<K, u128> for IntKmer<K, u128> {
    fn rev_comp(self) -> Self {
        let mut res = !self.to_int();
        res = (res >> 2 & 0x3333_3333_3333_3333_3333_3333_3333_3333)
            | (res & 0x3333_3333_3333_3333_3333_3333_3333_3333) << 2;
        res = (res >> 4 & 0x0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F)
            | (res & 0x0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F) << 4;
        res = (res >> 8 & 0x00FF_00FF_00FF_00FF_00FF_00FF_00FF_00FF)
            | (res & 0x00FF_00FF_00FF_00FF_00FF_00FF_00FF_00FF) << 8;
        res = (res >> 16 & 0x0000_FFFF_0000_FFFF_0000_FFFF_0000_FFFF)
            | (res & 0x0000_FFFF_0000_FFFF_0000_FFFF_0000_FFFF) << 16;
        res = (res >> 32 & 0x0000_0000_FFFF_FFFF_0000_0000_FFFF_FFFF)
            | (res & 0x0000_0000_FFFF_FFFF_0000_0000_FFFF_FFFF) << 32;
        res = (res >> 64 & 0x0000_0000_0000_0000_FFFF_FFFF_FFFF_FFFF)
            | (res & 0x0000_0000_0000_0000_FFFF_FFFF_FFFF_FFFF) << 64;
        Self::from_int(res >> (2 * (64 - K)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rc_8() {
        let kmer = IntKmer::<4, u8>::from_chars(b"ATCG");
        assert_eq!(kmer.rev_comp().to_chars(), *b"CGAT");
    }
    #[test]
    fn test_rc_16() {
        let kmer = IntKmer::<4, u16>::from_chars(b"ATCG");
        assert_eq!(kmer.rev_comp().to_chars(), *b"CGAT");
    }
    #[test]
    fn test_rc_32() {
        let kmer = IntKmer::<11, u32>::from_chars(b"CATAATCCAGC");
        assert_eq!(kmer.rev_comp().to_chars(), *b"GCTGGATTATG");
    }
    #[test]
    fn test_rc_64() {
        let kmer = IntKmer::<11, u64>::from_chars(b"CATAATCCAGC");
        assert_eq!(kmer.rev_comp().to_chars(), *b"GCTGGATTATG");
    }
    #[test]
    fn test_rc_128() {
        let kmer = IntKmer::<11, u128>::from_chars(b"CATAATCCAGC");
        assert_eq!(kmer.rev_comp().to_chars(), *b"GCTGGATTATG");
    }
}
