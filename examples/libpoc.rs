use argh::FromArgs;
use std::path::PathBuf;
use std::time::Instant;
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
    let start = Instant::now();

    let args: LibPoc = argh::from_env();

    // Read peptides
    ////////////////

    let pep_read_start = Instant::now();

    // Assume inputs are small and simply read into memory
    let pep_bytes = slurp_file(args.peptides_file);

    let pep_read_duration = pep_read_start.elapsed();
    println!("Read peptides file in {:.4} sec", pep_read_duration.as_secs_f64());

    // Build trie
    /////////////

    let trie_build_start = Instant::now();
    let mut trie = PeptideTrie::new();

    // Loop over peptides
    for pep in pep_bytes.split(|b| *b == b'\n') {
        if pep.len() == 0 {
            continue
        }

        trie.add(pep);
    }

    let trie_build_duration = trie_build_start.elapsed();
    println!("Added {} peptides to trie in {:.4} sec", trie.len(), trie_build_duration.as_secs_f64());

    // println!("STATS: {:#?}", trie.stats().expect("Trie was empty!"));
    //
    // println!("PEPTIDES: {:#?}", trie);

    // Annotate proteins
    ////////////////////

    let annotate_start = Instant::now();

    // let prots = annotate_fasta(args.fasta_file, trie);
    let prots = read_fasta(args.fasta_file);

    let mut num_entries : u64 = 0;
    // let mut total_edges : u64 = 0;
    for entry in prots.iter() {
        // println!("{}: {} peptides -- {:?}", entry.header(), entry.peptides().len(), entry.peptides());
        // println!("{}: {} peptides", entry.header(), entry.sequence());

        num_entries += 1;
        // total_edges += entry.peptides().len() as u64;
    }

    let annotate_duration = annotate_start.elapsed();
    // println!("Read and annotated {num_entries} entries with {total_edges} edges in {:.4} sec", annotate_duration.as_secs_f64());
    println!("Read and annotated {num_entries} entries in {:.4} sec", annotate_duration.as_secs_f64());
    println!("Total execution time: {:.4} sec", start.elapsed().as_secs_f64());
}