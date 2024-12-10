use argh::FromArgs;
use std::path::PathBuf;

use preppers::io::slurp_file;
use preppers::*;
use preppers::fasta::*;


/// Try building and using a trie of peptides
#[derive(FromArgs)]
struct LibPoc {
    /// path to peptides file
    #[argh(positional)]
    peptides_file: PathBuf,

    /// path to FASTA
    #[argh(positional)]
    fasta_file: PathBuf,
}


fn main() {
    let args: LibPoc = argh::from_env();

    // Assume inputs are small and simply read into memory
    let pep_bytes = slurp_file(args.peptides_file);

    let mut trie = PeptideTrie::new();

    // Loop over peptides
    for pep in pep_bytes.split(|b| *b == b'\n') {
        if pep.len() == 0 {
            continue
        }

        trie.add(pep);
    }

    println!("Added {} peptides to trie", trie.len());

    println!("STATS: {:#?}", trie.stats().expect("Trie was empty!"));

    let prots = annotate_fasta(args.fasta_file, trie);
    // let prots = read_fasta(args.fasta_file);

    let mut num_entries : u64 = 0;
    for entry in prots {
        println!("{}: {} peptides -- {:?}", entry.header(), entry.peptides().len(), entry.peptides());
        // println!("{}: {} peptides", entry.header(), entry.sequence());

        num_entries += 1;
    }
    println!("Read {num_entries} entries!")
}