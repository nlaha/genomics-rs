use std::{fs, process::exit};

use log::{error, info};
use serde_derive::Deserialize;


#[derive(Deserialize, Clone)]
#[serde(rename = "scores")]
pub struct Scores {
    pub s_match: i64,
    pub s_mismatch: i64,
    pub g: i64,
    pub h: i64,
}

#[derive(Deserialize)]
#[serde(rename = "fasta")]
pub struct FastaConfig {
    pub path: String,
}

#[derive(Deserialize)]
pub struct Config {
    pub scores: Scores,
    pub fasta: FastaConfig,
}

/// Read the config file and return a Config struct.
pub fn get_config(filepath: &str) -> Config
{
    // read config file contents to string
    let contents = match fs::read_to_string(filepath.to_string()) {
        Ok(c) => c,
        Err(_) => {
            error!("Could not read config file: {}", filepath);
            exit(1);
        }
    };

    let config: Config = match toml::from_str(&contents) {
        Ok(d) => d,
        Err(_) => {
            error!("Could not parse config file: {}", filepath);
            exit(1);
        }
    };

    // log config
    info!("Using the following values for scoring:");
    info!("Match: {}", config.scores.s_match);
    info!("Mismatch: {}", config.scores.s_mismatch);
    info!("Gap: {}", config.scores.g);
    info!("Opening Gap: {}", config.scores.h);

    info!("Using the following FASTA file: {}", config.fasta.path);

    return config;
}