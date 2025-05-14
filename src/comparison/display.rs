use colored::Colorize;

fn percentage(num: usize, total: usize) -> f64 {
    (num as f64 / total as f64) * 100.0
}

pub fn print_similarity_matrix(similarity_matrix: &ndarray::Array2<(usize, usize, usize, usize)>) {
    // print columns
    print!("  ");
    for i in 0..similarity_matrix.ncols() {
        print!("{} ", i);
    }
    println!();

    // display  table
    for i in 0..similarity_matrix.nrows() {
        // print sequence base
        print!("{} ", i);
        for j in 0..similarity_matrix.ncols() {
            let score = similarity_matrix[[i, j]];
            let max = std::cmp::max(score.1, score.2);
            let percentage_value = percentage(score.0, max);
            let color = VIRIDIS_COLORS[percentage_value as usize / 4];
            print!("{} ", "■".truecolor(color.0, color.1, color.2));
        }
        println!();
    }
}

/*
    This part of the code is generated using a GenAI tool. The prompt I used to generate this code is as follows:
        "Create a rust constant array with values (u8,u8,u8) representing rgb where the index maps to the corresponding
        viridis scale color. There should be a total of 26 entries"
    I bear full responsibility for this code, and I can explain its full behavior as needed.
    If I cannot explain this code satisfactorily, I risk losing points or potentially failing the assignment.
*/
pub const VIRIDIS_COLORS: [(u8, u8, u8); 26] = [
    (68, 1, 84),     // 0 - Dark purple
    (72, 22, 100),   // 1
    (71, 42, 113),   // 2
    (66, 63, 122),   // 3
    (59, 81, 128),   // 4
    (51, 99, 132),   // 5
    (43, 115, 134),  // 6
    (36, 131, 133),  // 7
    (31, 147, 129),  // 8
    (33, 163, 124),  // 9
    (42, 178, 116),  // 10
    (57, 192, 105),  // 11
    (74, 205, 93),   // 12
    (93, 217, 81),   // 13
    (114, 228, 69),  // 14
    (135, 238, 57),  // 15
    (157, 246, 47),  // 16
    (178, 253, 38),  // 17
    (199, 253, 33),  // 18
    (218, 251, 33),  // 19
    (234, 247, 34),  // 20
    (244, 241, 39),  // 21
    (249, 231, 47),  // 22
    (252, 220, 58),  // 23
    (253, 231, 37),  // 24 - yellow
    (253, 253, 253), // 25 - White
];
/* end of GenAI code */
