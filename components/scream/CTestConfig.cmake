
set(CTEST_SITE $ENV{CIME_MACHINE})

set(CTEST_PROJECT_NAME "SCREAM")
string(TIMESTAMP CURRTIME "%H:%M:%S" UTC)
set(CTEST_NIGHTLY_START_TIME "${CURRTIME} UTC")

set(CTEST_DROP_METHOD "https")
set (
  CTEST_CURL_OPTIONS
    "CURLOPT_SSL_VERIFYPEER_OFF"
    "CURLOPT_SSL_VERIFYHOST_OFF"
)
set(CTEST_DROP_SITE "my.cdash.org")
set(CTEST_DROP_LOCATION "/submit.php?project=E3SM")
set(CTEST_DROP_SITE_CDASH TRUE)
set(CTEST_COVERAGE_COMMAND "gcov")

set(CTEST_TEST_TIMEOUT 82800 CACHE STRING "")
set(DART_TESTING_TIMEOUT 82800 CACHE STRING "")

set(PARALLEL_LEVEL 1)
