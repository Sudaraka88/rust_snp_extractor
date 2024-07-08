use seq_io::{
    fasta::{write_head, write_seq, Reader},
    BaseRecord,
};

use std::{fs::File, io::Write, path::PathBuf};

use crate::{extract_sites::SnpPos, validate_inputs::Inputs};

#[derive(Debug, thiserror::Error)]
pub enum FastaError {
    #[error("Please provide a valid fasta file as input")]
    ErrorOpeningFasta,
    #[error("Write folder exists, plesae provide a new folder to save outputs")]
    ErrorReadingFasta,
    #[error("Sequence length mismatch: Previous sequence(s) had length: {0}, but sequence: {1} has a length of {2}")]
    SequenceMismatch(usize, String, usize),
    #[error("Cannot create write file: {0}")]
    ErrorCreatingOutFasta(String),
    #[error("Cannot write head: {0}")]
    ErrorWritingHead(String),
    #[error("Cannot write sequence: {0}")]
    ErrorWritingSequence(String),
    #[error("Cannot write position to file: {0}")]
    ErrorWritingPos(usize),
}

#[derive(Debug)]
pub struct FastaFile {
    fasta_sites: Vec<u8>,
    seqlen: usize,
    num_seqs: u16,
    headers: Vec<Vec<u8>>,
}

impl FastaFile {
    pub fn read_fasta(filepath: PathBuf) -> Result<FastaFile, FastaError> {
        // let file = File::open(filepath).expect("Error opening fasta file");
        let file = File::open(&filepath).map_err(|_| FastaError::ErrorOpeningFasta)?;

        let mut reader = Reader::new(file);
        let mut fasta_sites: Vec<u8> = Vec::new();
        let mut first_record = true;
        let mut seqlen = 0;
        let mut num_seqs: u16 = 0;
        let mut headers = Vec::new();

        while let Some(result) = reader.next() {
            // let record = result.expect("Reading error");
            let record = result.map_err(|_| FastaError::ErrorReadingFasta)?;

            if first_record {
                seqlen = record.seq_lines().fold(0, |_l, seq| 1 + seq.len());
                first_record = false;
            } else {
                let seqlen_t = record.seq_lines().fold(0, |_l, seq| 1 + seq.len());
                if seqlen != seqlen_t {
                    return Err(FastaError::SequenceMismatch(
                        seqlen - 1,
                        record.id().unwrap().to_string(),
                        seqlen_t - 1,
                    ));
                }
            }

            // Keep appending to fasta_sites (converting the fasta to a vector!)
            fasta_sites.extend_from_slice(record.seq());

            // Collect the fasta headers to a vector (u8)
            headers.push(record.head().to_vec());
            // println!("<<{:?}>>", x)

            num_seqs += 1;
        }
        Ok(FastaFile {
            fasta_sites,
            seqlen,
            num_seqs,
            headers,
        })
        // println!("<<{:?}>>", fasta_sites)
    }

    pub fn fasta_sites(&self) -> &Vec<u8> {
        &self.fasta_sites
    }

    pub fn seqlen(&self) -> usize {
        self.seqlen
    }

    pub fn num_seqs(&self) -> u16 {
        self.num_seqs
    }

    pub fn headers(&self) -> &Vec<Vec<u8>> {
        &self.headers
    }

    pub fn write_fasta(&self, snp_pos: &SnpPos, inputs: &Inputs) -> Result<(), FastaError> {
        // Let's write the positions first
        let out_pos_path = inputs.write_folder().join("snps.pos");
        let mut file_pos = File::create(&out_pos_path).map_err(|_| {
            FastaError::ErrorCreatingOutFasta(out_pos_path.to_string_lossy().to_string())
        })?;

        for &num in snp_pos.snp_pos() {
            writeln!(file_pos, "{}", num + 1).map_err(|_| FastaError::ErrorWritingPos(num))?;
            // offset by 1 for 1-based positions
        }

        // Now write the fasta file
        let out_fasta_path = inputs.write_folder().join("snps.fa");
        let file = File::create(&out_fasta_path).map_err(|_| {
            FastaError::ErrorCreatingOutFasta(out_fasta_path.to_string_lossy().to_string())
        })?;

        for i in 0..self.headers.len() {
            let head = self.headers()[i].as_slice();
            write_head(&file, head).map_err(|_| {
                FastaError::ErrorWritingHead(String::from_utf8(head.to_vec()).unwrap())
            })?; // header

            // for i in 0..num_seqs {
            let index: Vec<usize> = snp_pos
                .snp_pos()
                .iter()
                .map(|&x| x + i * self.num_seqs as usize)
                .collect(); // position indexes

            let vals: Vec<u8> = index
                .iter()
                .filter_map(|&idx| self.fasta_sites.get(idx).copied())
                .collect(); // position values (we have to copy them because &[u8] needs to be contiguous)

            write_seq(&file, &vals).map_err(|_| {
                FastaError::ErrorWritingSequence(String::from_utf8(head.to_vec()).unwrap())
            })?;
        }
        Ok(())
    }
}
