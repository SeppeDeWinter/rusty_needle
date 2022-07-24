
use std::collections::HashMap;

struct Alignment {
    seq_a: String,
    seq_b: String,
    indel_score: i32,
    match_score: i32,
    mismatch_score: i32,
    alignment_grid: Vec<Vec<i32>>
}

fn init_alignment_grid(
    seq_a: &String, 
    seq_b: &String
) -> Vec<Vec<i32>> {
    let mut alignment_grid: Vec<Vec<i32>> = vec![
        vec![0; seq_a.chars().count() + 1]; seq_b.chars().count() + 1];
    
    for i in 0..seq_a.chars().count() + 1 {
        alignment_grid[0][i] = -(i as i32);
    }

    for j in 0..seq_b.chars().count() + 1 {
        alignment_grid[j][0] = -(j as i32);
    }
    alignment_grid
}

fn fill_alignment_grid(
    mut aligment: Alignment,
) -> Alignment{
    for i in 1..aligment.seq_a.chars().count() + 1 {
        for j in 1..aligment.seq_b.chars().count()  + 1 {
            let current_seq_a_char = aligment.seq_a.chars().nth(i - 1);
            let current_seq_b_char = aligment.seq_b.chars().nth(j - 1);
            let diagonal_score = if current_seq_a_char == current_seq_b_char {
                aligment.match_score
            } else {
                aligment.mismatch_score
            };
            let possible_scores = vec![
                aligment.alignment_grid[j    ][i - 1] + aligment.indel_score,
                aligment.alignment_grid[j - 1][i    ] + aligment.indel_score,
                aligment.alignment_grid[j - 1][i - 1] + diagonal_score,           
            ];
            aligment.alignment_grid[j][i] = *possible_scores.iter().max().unwrap();
        }
    }
    aligment
}

fn traverse_alignment_grid(
    alignment: &Alignment,
    current_pos: (usize, usize),
    delta_pos: (usize, usize),
    aligned_seq_a: &String,
    aligned_seq_b: &String
) -> Vec<(String, String)> {

    let aligned_seq_b = if delta_pos.0 == 0 {
        "-".to_string() + &aligned_seq_b
    } else {
        alignment.seq_b.chars().nth(current_pos.0).unwrap().to_string() + &aligned_seq_b
    };

    let aligned_seq_a = if delta_pos.1 == 0 {
        "-".to_string() + &aligned_seq_a
    } else {
        alignment.seq_a.chars().nth(current_pos.1).unwrap().to_string() + &aligned_seq_a
    };

    let current_score = alignment.alignment_grid[current_pos.0][current_pos.1];
    
    let current_seq_a_char = alignment.seq_a.chars().nth(current_pos.1 - 1);
    let current_seq_b_char = alignment.seq_b.chars().nth(current_pos.0 - 1);
    let diagonal_score = if current_seq_a_char == current_seq_b_char {
        alignment.match_score
    } else {
        alignment.mismatch_score
    };

    let mut routes = HashMap::new();
    routes.insert(
        (current_pos.0, current_pos.1 - 1), 
        alignment.alignment_grid[current_pos.0][current_pos.1 - 1] + alignment.indel_score);
    routes.insert(
        (current_pos.0 - 1, current_pos.1), 
        alignment.alignment_grid[current_pos.0 - 1][current_pos.1] + alignment.indel_score);
    routes.insert(
        (current_pos.0 - 1, current_pos.1 - 1), 
        alignment.alignment_grid[current_pos.0 - 1][current_pos.1 - 1] + diagonal_score);
    routes.retain(|_, &mut v| v == current_score);
    let routes = routes.keys();
    
    let mut return_val = Vec::new();
    for route in routes {
        if route != &(0, 0) {
            let delta_pos = (current_pos.0 - route.0, current_pos.1 - route.1);
            let current_pos = (current_pos.0 - delta_pos.0, current_pos.1 - delta_pos.1);
            return_val.append(&mut traverse_alignment_grid(
                &alignment,
                current_pos,
                delta_pos,
                &aligned_seq_a,
                &aligned_seq_b))

        }
        else {
            let mut aligned_seq_b = alignment.seq_b.chars().nth(0).unwrap().to_string() + &aligned_seq_b;
            let mut aligned_seq_a = alignment.seq_a.chars().nth(0).unwrap().to_string() + &aligned_seq_a;
            //remove dash from end of the sequence
            aligned_seq_b.pop().unwrap();
            aligned_seq_a.pop().unwrap();
            return vec![(aligned_seq_a, aligned_seq_b)]
        }
    }
    return return_val
}

pub fn align(seq_a: String, seq_b: String) {
    let alignment_grid = init_alignment_grid(&seq_a, &seq_b);
    let mut alignment = Alignment{
        seq_a: seq_a, 
        seq_b: seq_b,
        indel_score: -1,
        mismatch_score: -1,
        match_score: 1,
        alignment_grid: alignment_grid  
    };
    alignment = fill_alignment_grid(alignment);
    let current_pos =  (alignment.seq_b.chars().count(), alignment.seq_a.chars().count());
    let alignments = traverse_alignment_grid(
        &alignment,
        current_pos,
        (0, 0),
        &"".to_string(),
        &"".to_string()
    );
    for (i, (aligned_seq_a, aligned_seq_b)) in alignments.iter().enumerate() {
        println!("Alignment {}", i + 1);
        println!("seq_a: {}", aligned_seq_a);
        println!("seq_b: {}", aligned_seq_b);
        println!("")
    }
}