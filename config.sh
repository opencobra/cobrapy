# Define custom utilities
# Test for OSX with [ -n "$IS_OSX" ]

function pre_build {
    # Any stuff that you need to do before you start building the wheels
    # Runs in the root directory of this repository.
    if [ -n "$IS_OSX" ]; then
        export CC=clang
        export CXX=clang++
		export CFLAGS="-fPIC -O3 -arch i386 -arch x86_64 -g -DNDEBUG -mmacosx-version-min=10.6"
	else
		yum install -y libxslt libxml2 libxml2-devel libxslt-devel
	fi
	curl -O http://ftp.gnu.org/gnu/glpk/glpk-4.60.tar.gz
	tar xzf glpk-4.60.tar.gz
	(cd glpk-4.60 \
			&& ./configure \
			&& make \
			&& make install)
	pip install cython
	cython -a cobra/solvers/cglpk.pyx
	export PATH="/usr/local/bin:$PATH"
}

function build_wheel {
    # Set default building method to pip
    build_bdist_wheel $@
}

function run_tests_in_repo {
    # Run tests from within source repo
	# trick is to run the /installed/ package 
    # coverage run --source=cobra setup.py test
	if [ -n "$IS_OSX" ]; then
		echo -e " testing for mac.. "
	else
		wget --no-check-certificate \
			 https://opencobra.github.io/pypi_cobrapy_travis/esolver.gz
		echo -e " travis randomly stopping? 1"
		gzip -f -d esolver.gz
		chmod +x esolver
		export PATH=$PATH:$PWD
		echo -e " travis randomly stopping? 2"
		which pkg-config
		pip install matplotlib
		echo -e " travis randomly stopping? 3"
	fi
	echo -e " travis randomly stopping? 4"
	mkdir -p $HOME/.config/matplotlib
	echo 'backend: Agg' >> $HOME/.config/matplotlib/matplotlibrc
	echo -e " travis randomly stopping? 5"
	python -c "import cobra.test; cobra.test.test_all()"
	echo $?
}

function run_tests {
    # Runs tests on installed distribution from an empty directory
    run_tests_in_repo
}
