
pub(crate) const MAXIMUM_K_SIZE: usize = u32::max_value() as usize;
#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("K size {ksize} is out of range for the given sequence size {seq_size}")]
    KSizeOutOfRange { ksize: usize, seq_size: usize },
    #[error("K size {0} cannot exceed the size of a u32 {MAXIMUM_K_SIZE}")]
    KSizeTooBig(usize),
}

pub type Result<T> = std::result::Result<T, Error>;
