use std::fs::OpenOptions;
use std::io::Read;
use std::path::PathBuf;

pub fn slurp_file(input_file: PathBuf) -> Vec<u8> {
    let mut input_file = OpenOptions::new()
        .read(true)
        .open(input_file)
        .expect("unable to open text file");

    let metadata = input_file
        .metadata()
        .expect("file metadata should be available");
    let file_len = metadata.len() as usize;
    let mut contents = Vec::with_capacity(file_len);

    let bytes_read = input_file
        .read_to_end(&mut contents)
        .expect("should be able to read entire file");

    if file_len > 0 {
        assert_eq!(
            bytes_read, file_len,
            "should have read entire file into memory"
        );
    }

    contents
}
