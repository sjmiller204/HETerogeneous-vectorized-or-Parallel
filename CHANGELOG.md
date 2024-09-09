# Changelog

## [Unreleased] - TBD
### Changed
- Moved fortran files to src/Core and src/Test
- Changed hetp_main.F90 to print test output to file in same directory as executable

## Added
- Added CMake files to enable cmake build
- Added CHANGELOG.md to record changes per version
- Added .gitignore to avoid committing temporary files
- Added REAMDE instructions to run standalone test and connect to external models
- Added github action to test build and display badge in the README
- Added option to write test inputs and outputs to terminal (commented out by default)

### Fixed
- Fixed initialization bugs found using compiler debug flags
