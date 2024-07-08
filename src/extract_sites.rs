use crate::fasta_io::FastaFile;

#[derive(Debug, thiserror::Error)]
pub enum ExtractionError {
    #[error("Your MSA has no SNPs")]
    EmptySNPVector,
}

pub struct SnpPos {
    snp_pos: Vec<usize> // Add functionality to extract allele frequencies
}

impl SnpPos{
    pub fn get_snp_positions(fastafile: &FastaFile, maf: f32) -> Result<SnpPos, ExtractionError> {
        // let end = (seq_len - 1)*(num_seqs as usize);
        // let snp: Vec<u8> = Vec::new();
        // let mut k = 0;
        let mut snp_pos = Vec::new();
        let min_alleles: i32 = (maf * (fastafile.num_seqs() as f32)).ceil() as i32;
        
        // println!("ma:{}", min_alleles);

        // This can probably be parallelised    
        // This is essentially a linearised structure scanned using pseudo row/column indices
        for i in 0..(fastafile.seqlen()-1) { // 0 indexing
            // println!("i:{}", i);
            // row index
            let mut num_alleles: i32 = 0; // How many alleles in this SNP?
    
            for j in 1..(fastafile.num_seqs() as usize) {
                // column index
                // println!("idx:{}", j * (fastafile.seqlen()-1) + i);
                if fastafile.fasta_sites()[i] != fastafile.fasta_sites()[j * (fastafile.seqlen()-1) + i] {
                    num_alleles += 1;
                    // println!("na:{}", num_alleles);
    
                    if num_alleles >= min_alleles {
                        snp_pos.push(i);
                        break;
                    }
                }
            }
        }
        if snp_pos.len() == 0 {
            return Err(ExtractionError::EmptySNPVector)
        } else {
            Ok(SnpPos{ snp_pos })
        }
    }
    
    pub fn snp_pos(&self) -> &[usize] {
        &self.snp_pos
    } 

}
