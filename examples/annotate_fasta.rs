use argh::FromArgs;
use preppers::fasta::*;
use preppers::io::slurp_file;
use preppers::*;
use std::path::PathBuf;
use std::time::Instant;


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
    let mut n = 0;
    for pep in pep_bytes.split(|b| *b == b'\n') {
        if pep.len() == 0 {
            continue
        }

        trie.insert(pep);

        n += 1;
    }

    let trie_build_duration = trie_build_start.elapsed();
    println!("Read {} sequences and built trie for {} peptides in {:.4} sec", n, trie.len(), trie_build_duration.as_secs_f64());

    // println!("STATS: {:#?}", trie.stats().expect("Trie was empty!"));
    collect_and_output_stats(&trie);

    // println!("PEPTIDES: {:#?}", trie);

    // Annotate proteins
    ////////////////////

    let annotate_start = Instant::now();

    let fasta_read_start = Instant::now();

    let fasta = read_fasta(args.fasta_file);

    let fasta_read_duration = fasta_read_start.elapsed();
    println!("Read FASTA in {:.4} sec", fasta_read_duration.as_secs_f64());

    // let prots = fasta.iter();
    let prots = annotate_fasta(&fasta, trie).unwrap();

    let mut num_entries : u64 = 0;
    let mut total_edges : u64 = 0;

    for entry in prots {
        // println!("{}: {} peptides -- {:?}", entry.header(), entry.peptides().len(), entry.peptides());
        // println!("{}: {} peptides", entry.header(), entry.sequence());

        num_entries += 1;
        total_edges += entry.peptides().len() as u64;
    }

    let annotate_duration = annotate_start.elapsed();
    // println!("Read and annotated {num_entries} entries in {:.4} sec", annotate_duration.as_secs_f64());
    println!("Parsed and annotated {num_entries} entries with {total_edges} edges in {:.4} sec", annotate_duration.as_secs_f64());
    println!("Total execution time: {:.4} sec", start.elapsed().as_secs_f64());
}

fn collect_and_output_stats(peptides: &PeptideTrie) -> Option<()> {
    let stats = peptides.stats()?;

    println!("STATS: {stats}");

    let overhead_bytes_per_key_byte =
        (stats.tree.mem_usage as f64) / (stats.leaf.sum_key_bytes as f64);

    println!("{overhead_bytes_per_key_byte} bytes of overhead, per byte of key stored in tree");

    Some(())
}