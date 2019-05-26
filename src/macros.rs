/// Takes a collection of collections and pads them with the given element such that all collections are of same length.
/// Returns the new length. Mutates the collections.
#[macro_export]
macro_rules! pad_collection {
    ( $coll:expr, $pad:expr ) => {
        {
            if $coll.len() > 0 {
                // Find maximum length of a collection
                let mut max_length = $coll[0].len();
                for i in 1..$coll.len() {
                    if max_length < $coll[i].len() {
                        max_length = $coll[i].len();
                    }
                }
                for i in 0..$coll.len() {
                    let l = $coll[i].len();
                    if l < max_length {
                        // append might not unroll and hence perform slower
                        //$coll[i].append(&mut vec![$pad; max_length - l]);
                        for _ in 0..max_length - l {
                            $coll[i].push($pad);
                        }
                    }
                }

                max_length
            } else {
                0
            }
        }
    };
}
