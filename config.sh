# Define custom utilities
# Test for OSX with [ -n "$IS_OSX" ]

function pre_build {
    # Any stuff that you need to do before you start building the wheels
    # Runs in the root directory of this repository.
	wget --no-check-certificate https://opencobra.github.io/pypi_cobrapy_travis/esolver.gz
	gzip -f -d esolver.gz		# download seems to happen twice? so force
	chmod +x esolver
	export PATH=$PATH:$PWD
	mkdir -p $HOME/.config/matplotlib
	echo 'backend: Agg' >> $HOME/.config/matplotlib/matplotlibrc
    if [ -n "$IS_OSX" ]; then
        export CC=clang
        export CXX=clang++
		export CFLAGS="-fPIC -O3 -arch i386 -arch x86_64 -g -DNDEBUG -mmacosx-version-min=10.6"
		brew tap homebrew/science
        brew update
        brew install glpk
	else
		curl -O http://ftp.gnu.org/gnu/glpk/glpk-4.60.tar.gz
		tar xzf glpk-4.60.tar.gz
		(cd glpk-4.60 \
				&& ./configure --prefix=$BUILD_PREFIX \
				&& make \
				&& make install)
		cp $BUILD_PREFIX/include/glpk.h .
		cp $BUILD_PREFIX/lib/libglpk.a .
		yum install -y libxslt libxml2 libxml2-devel libxslt-devel
	fi
	pip install cython
	cython -a cobra/solvers/cglpk.pyx
}

function build_wheel {
    # Set default building method to pip
    build_bdist_wheel $@
}

function run_tests_in_repo {
    # Run tests from within source repo
	ls -la
	pwd
	echo $BUILD_PREFIX
	ls -la $BUILD_PREFIX
    coverage run --source=cobra setup.py test
}

function run_tests {
    # Runs tests on installed distribution from an empty directory
    (cd .. && run_tests_in_repo)
}
