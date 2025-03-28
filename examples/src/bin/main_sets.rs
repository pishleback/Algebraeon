use algebraeon::sets::combinatorics::*;

fn main() {
    let mut c = LexicographicSubsetsWithRemovals::new(7, 3);
    for _ in 0..19 {
        let x = c.next().unwrap();
        println!("{:?}", x);
    }

    println!("rm 4");
    c.exclude(4);
    println!("rm 5");
    c.exclude(5);

    for _ in 0..2 {
        let x = c.next().unwrap();
        println!("{:?}", x);
    }

    assert_eq!(c.next(), None);
}
