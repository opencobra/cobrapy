# Define custom utilities
# Test for OSX with [ -n "$IS_OSX" ]

function pre_build {
    # Any stuff that you need to do before you start building the wheels
    # Runs in the root directory of this repository.
	wget https://opencobra.github.io/pypi_cobrapy_travis/esolver.gz
	gzip -d esolver.gz
	chmod +x esolver
	export PATH=$PATH:$PWD
	mkdir -p $HOME/.config/matplotlib
	echo 'backend: Agg' >> $HOME/.config/matplotlib/matplotlibrc
    curl -O http://ftp.gnu.org/gnu/glpk/glpk-4.60.tar.gz
    tar xzf glpk-4.60.tar.gz
    if [ -n "$IS_OSX" ]; then
        export CC=clang
        export CXX=clang++
		export CFLAGS="-fPIC -O3 -arch x86_64 -g -DNDEBUG -mmacosx-version=10.6"
	fi
	(cd glpk-4.60 \
			&& ./configure --prefix=$BUILD_PREFIX \
			&& make \
			&& make install)
}

function build_wheel {
    # Set default building method to pip
    build_bdist_wheel $@
}

function run_tests_in_repo {
    # Run tests from within source repo
    coverage run --source=cobra setup.py test
}

function run_tests {
    # Runs tests on installed distribution from an empty directory
    run_tests_in_repo
}
