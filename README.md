# genomics_rs

This is the repository for my CPT_S 471 Computational Genomics homework. This application, written in Rust, aligns FASTA sequences via local and global alignment with the affine gap penalty function.

To run the program, configure parameters and input file in `config.toml`, then use `cargo run --release`. While running it without the `--release` flag is possible, performance will drop substantially.

There is visualization support for smaller sequences:
![image](https://github.com/user-attachments/assets/57e57f96-cc44-43c0-821f-5e6b3839a84a)
