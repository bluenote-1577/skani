#!/bin/bash

# Define the expected version
EXPECTED_VERSION="0.3.0"

# Function to extract the version from Cargo.toml
get_cargo_version() {
      grep -m1 "^version" Cargo.toml | sed -E 's/version\s*=\s*"([^"]+)"/\1/'
}

get_main_version() {
      grep 'pub const VERSION' src/params.rs | sed -E 's/pub const VERSION\s*:\s*&str\s*=\s*"([^"]+)";/\1/'
}

# Check if Cargo.toml has the expected version
CARGO_VERSION=$(get_cargo_version)
if [ "$CARGO_VERSION" == "$EXPECTED_VERSION" ]; then
  echo "Cargo.toml has the correct version: $CARGO_VERSION."
else
  echo "Error: Cargo.toml version ($CARGO_VERSION) does not match expected version ($EXPECTED_VERSION)."
  exit 1
fi

# Check if main.rs has the expected version
MAIN_VERSION=$(get_main_version)
if [ "$MAIN_VERSION" == "$EXPECTED_VERSION" ]; then
  echo "main.rs has the correct version: $MAIN_VERSION."
else
  echo "Error: main.rs version ($MAIN_VERSION) does not match expected version ($EXPECTED_VERSION)."
  exit 1
fi

# Run the cargo test command
echo "Running tests..."
cargo test -j 1 -- --show-output > test_results_versions/$EXPECTED_VERSION

echo "Test results have been saved to test_results_versions/$EXPECTED_VERSION"
