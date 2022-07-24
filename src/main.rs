fn main() {
    let seq_a = String::from("ATTCTATCAAGGTAGCTTGATTCCAGGTCGATTCC");
    let seq_b = String::from("ATTCTAGCAAGGTAGCTTGATTCCAGTCGATTCC");

    rusty_needle::align(seq_a, seq_b)
}
