docker build -t cobrapy_builder .
docker run --rm -v `pwd`:/io cobrapy_builder /io/build_cobrapy.sh
