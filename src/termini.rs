use fancy_regex::Regex;
use crate::fasta::{PeptideHit, PreppedFastaEntry};

pub fn filter_entry_termini<'a, T: Iterator<Item=PreppedFastaEntry<'a>>>(iter: T, enzyme_patt: Regex, n_req_termini: u8) -> impl Iterator<Item=PreppedFastaEntry<'a>> {
    iter.map(move |entry| filter_match_termini(entry, &enzyme_patt, n_req_termini))
}

fn filter_match_termini<'a>(mut entry: PreppedFastaEntry<'a>, enzyme_patt: &Regex, n_req_termini: u8) -> PreppedFastaEntry<'a> {
    if n_req_termini > 0 {
        let seq = entry.sequence.as_slice();

        entry
            .peptide_indices
            .retain(|hit| {
                has_required_termini(hit, seq, enzyme_patt, n_req_termini)
            });
    }

    entry
}

pub fn has_required_termini(peptide_hit: &PeptideHit, seq: &[u8], enzyme_patt: &Regex, n_req_termini: u8) -> bool {
    match n_req_termini {
        0 => true,
        _ => num_termini(peptide_hit, seq, enzyme_patt) >= n_req_termini
    }
}

pub fn num_termini(peptide_hit: &PeptideHit, seq: &[u8], enzyme_patt: &Regex) -> u8 {
    let (_, start, stop) = *peptide_hit;

    let seq_str = match std::str::from_utf8(seq) {
        Ok(s) => s,
        Err(_) => return 0,
    };

    let boundary_match = |idx: usize| -> bool {
        if idx == 0 || idx == seq.len() {
            return true;
        }

        // If the given index isn't as character boundary we can't match it
        if !seq_str.is_char_boundary(idx) {
            return false;
        }

        match enzyme_patt.find_from_pos(seq_str, idx) {
            Ok(Some(m)) => m.start() == idx,
            _ => false,
        }
    };

    // Note: we check for match at `stop + 1` because the zero-width match should occur _after_ the
    // last character of the peptide.
    boundary_match(start) as u8 + boundary_match(stop + 1) as u8
}

#[cfg(test)]
mod tests {
    use fancy_regex::Regex;
    use crate::termini::{has_required_termini, num_termini};

    #[test]
    fn test_num_termini_counts_matches() {
        let seq = b"APEPTIDEKANOTHER";
        let enzyme_patt = Regex::new(r"(?<=[KR])").unwrap();

        assert_eq!(num_termini(&(999, 9, 15), seq, &enzyme_patt), 2);
        assert_eq!(num_termini(&(999, 10, 15), seq, &enzyme_patt), 1);
    }

    #[test]
    fn test_has_required_termini_matches_exact_positions() {
        let seq = b"ABCDE";

        // Match peptides starting with C or ending with E
        let enzyme_patt = Regex::new(r"(?=C)|(?<=E)").unwrap();

        assert!(has_required_termini(&(999, 2, 4), seq, &enzyme_patt, 2));
        assert!(!has_required_termini(&(999, 1, 4), seq, &enzyme_patt, 2));
    }

    #[test]
    fn test_has_required_termini_accepts_sequence_boundaries() {
        let seq = b"ABCDE";

        // Match peptides ending before an E
        let enzyme_patt = Regex::new(r"(?=E)").unwrap();

        assert!(has_required_termini(&(999, 0, 3), seq, &enzyme_patt, 2));
    }

    #[test]
    fn test_num_termini_counts_termini() {
        let seq = b"ABCDE";

        // Match peptides ending before a Z
        let enzyme_patt = Regex::new(r"(?=Z)").unwrap();

        // Expect one match, from the peptide starting at the N-terminus
        assert_eq!(num_termini(&(999, 0, 4), seq, &enzyme_patt), 2);
    }

    #[test]
    fn test_num_termini_counts_n_term() {
        let seq = b"ABCDE";

        // Match peptides ending before a Z
        let enzyme_patt = Regex::new(r"(?=Z)").unwrap();

        // Expect one match, from the peptide starting at the N-terminus
        assert_eq!(num_termini(&(999, 0, 3), seq, &enzyme_patt), 1);
    }

    #[test]
    fn test_num_termini_counts_c_term() {
        let seq = b"ABCDE";

        // Match peptides ending before a Z
        let enzyme_patt = Regex::new(r"(?=Z)").unwrap();

        // Expect one match, from the peptide ending at the C-terminus
        assert_eq!(num_termini(&(999, 1, 4), seq, &enzyme_patt), 1);
    }
}
