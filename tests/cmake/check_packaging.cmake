# Checks that release metadata files exist and that the changelog is current.
# Inputs: SOURCE_DIR, PROJECT_VERSION
foreach(f CITATION.cff CHANGELOG.md)
  if(NOT EXISTS ${SOURCE_DIR}/${f})
    message(FATAL_ERROR "missing ${f}")
  endif()
endforeach()

file(READ ${SOURCE_DIR}/CHANGELOG.md changelog)
string(FIND "${changelog}" "${PROJECT_VERSION}" has_version)
string(FIND "${changelog}" "Unreleased" has_unreleased)
if(has_version EQUAL -1 AND has_unreleased EQUAL -1)
  message(FATAL_ERROR "CHANGELOG.md mentions neither ${PROJECT_VERSION} nor 'Unreleased'")
endif()

file(READ ${SOURCE_DIR}/CITATION.cff cff)
foreach(doi 10.1093/bioinformatics/btad019 10.1093/bioinformatics/btq614)
  string(FIND "${cff}" "${doi}" pos)
  if(pos EQUAL -1)
    message(FATAL_ERROR "CITATION.cff does not contain doi ${doi}")
  endif()
endforeach()
