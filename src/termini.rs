use fancy_regex::Regex;
use crate::fasta::{PeptideHit, PreppedFastaEntry};

pub fn filter_entry_termini<'a, T: Iterator<Item=PreppedFastaEntry<'a>>>(iter: T, enzyme_patt: Regex, n_req_termini: u8) -> impl Iterator<Item=PreppedFastaEntry<'a>> {
    iter.map(move |entry| filter_match_termini(entry, &enzyme_patt, n_req_termini))
}

fn filter_match_termini<'a>(mut entry: PreppedFastaEntry<'a>, enzyme_patt: &Regex, n_req_termini: u8) -> PreppedFastaEntry<'a> {
    let seq = entry.sequence.as_slice();

    entry
        .peptide_indices
        .retain(|hit| {
            has_required_termini(hit, seq, enzyme_patt, n_req_termini)
        });

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
        if idx == 0 || idx == seq.len() - 1 {
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

    boundary_match(start) as u8 + boundary_match(stop) as u8
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
        let enzyme_patt = Regex::new(r"(?=C)|(?=E)").unwrap();

        assert!(has_required_termini(&(999, 2, 4), seq, &enzyme_patt, 2));
        assert!(!has_required_termini(&(999, 1, 4), seq, &enzyme_patt, 2));
    }

    #[test]
    fn test_has_required_termini_accepts_sequence_boundaries() {
        let seq = b"ABCDE";
        let enzyme_patt = Regex::new(r"(?=E)").unwrap();

        assert!(has_required_termini(&(999, 0, 4), seq, &enzyme_patt, 2));
    }
}
