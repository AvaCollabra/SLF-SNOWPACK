
# add the tests
ADD_TEST(basics.smoke run_basics.sh)
SET_TESTS_PROPERTIES(basics.smoke
                     PROPERTIES LABELS smoke
                     FAIL_REGULAR_EXPRESSION "error|differ|fail")
