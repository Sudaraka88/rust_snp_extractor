use std::{fs, path::PathBuf};

#[derive(Debug, PartialEq, Clone)]
pub struct Inputs {
    read_path: PathBuf,
    write_folder: PathBuf,
    maf: f32,
}

#[derive(Debug, thiserror::Error)]
pub enum InputError {
    #[error("Please provide a valid fasta file as input")]
    WrongFasta,
    #[error("Write folder exists, plesae provide a new folder to save outputs")]
    WriteFolderExists,
    #[error("Minor allele frequency should be between (0,0.5)")]
    WrongMAF,
    #[error("Folder cannot be created at the provided path: {0}")]
    FolderCreationError(String),
}

impl Inputs {
    pub fn parse_inputs(
        read_path: String,
        write_folder: String,
        maf: f32,
    ) -> Result<Self, InputError> {
        let read_path = PathBuf::from(&read_path);
        if !read_path.exists() {
            return Err(InputError::WrongFasta);
        }

        let write_folder = PathBuf::from(&write_folder);
        if write_folder.exists() {
            return Err(InputError::WriteFolderExists);
        }

        fs::create_dir(&write_folder).map_err(|err| {
            InputError::FolderCreationError(format!(
                "Folder cannot be created at the provided path: {}",
                err
            ))
        })?;

        if maf < 0.0 || maf > 0.5 {
            return Err(InputError::WrongMAF);
        }

        Ok(Inputs {
            read_path,
            write_folder,
            maf,
        })
    }
    
    pub fn read_path(&self) -> &PathBuf {
        &self.read_path
    }
    
    pub fn write_folder(&self) -> &PathBuf {
        &self.write_folder
    }
    
    pub fn maf(&self) -> f32 {
        self.maf
    }
}
