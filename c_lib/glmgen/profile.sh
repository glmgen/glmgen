valgrind --tool=callgrind --callgrind-out-file=profile_raw --compress-strings=no --zero-before=tf_admm ./examples/bin/test_admm

callgrind_annotate --threshold=20 --context=1 --inclusive=yes --auto=yes profile_raw > profile.out
