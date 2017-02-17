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
    curl -O http://ftp.gnu.org/gnu/glpk/glpk-4.61.tar.gz
    tar xzf glpk-4.61.tar.gz
    (cd glpk-4.61 \
            && ./configure --disable-reentrant \
            && make \
            && make install)
    pip install cython
    cython -a cobra/solvers/cglpk.pyx
    export PATH="$PATH:/usr/local/bin"
}

function build_wheel {
    # Set default building method to pip
    build_bdist_wheel $@
    # since swiglpk doesn't have wheels, we currently must keep glpk
    # installed for testing
	# (cd glpk-4.60 && make uninstall)
}

function run_tests_in_repo {
    # Run tests from within source repo
	# trick is to run the /installed/ package 
    # coverage run --source=cobra setup.py test
	if [ -n "$IS_OSX" ]; then
		echo -e " testing for mac... "
	else
		wget --no-check-certificate \
			 https://opencobra.github.io/pypi_cobrapy_travis/esolver.gz
		gzip -f -d esolver.gz
		chmod +x esolver
		export PATH=$PWD:$PATH
		# which pkg-config
		# pip install matplotlib
	fi
	mkdir -p $HOME/.config/matplotlib
	echo 'backend: Agg' >> $HOME/.config/matplotlib/matplotlibrc
	COVERAGEXML=`python -c "import os,sys; print(os.path.realpath('coverage.xml'))"`
	COVERAGERC=`python -c "import os,sys; print(os.path.realpath('../.coveragerc'))"`
	(pytest --pyargs -v -rsx --cov=cobra --cov-report=xml:${COVERAGEXML} \
			--cov-config=${COVERAGERC} --benchmark-skip cobra &&
			mv ${COVERAGEXML} ..)
}

function run_tests {
    # Runs tests on installed distribution from an empty directory
    run_tests_in_repo
}
