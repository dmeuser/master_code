FILE(REMOVE_RECURSE
  "CMakeFiles/lt"
  "CMakeFiles/lt-complete"
  "lt-prefix/src/lt-stamp/lt-install"
  "lt-prefix/src/lt-stamp/lt-mkdir"
  "lt-prefix/src/lt-stamp/lt-download"
  "lt-prefix/src/lt-stamp/lt-update"
  "lt-prefix/src/lt-stamp/lt-patch"
  "lt-prefix/src/lt-stamp/lt-configure"
  "lt-prefix/src/lt-stamp/lt-build"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/lt.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
