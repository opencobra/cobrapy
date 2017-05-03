FROM quay.io/pypa/manylinux1_x86_64

ENV GLPK_VER="4.61"
RUN wget http://ftp.gnu.org/gnu/glpk/glpk-${GLPK_VER}.tar.gz -O - | tar xz
WORKDIR glpk-${GLPK_VER}
RUN ./configure && make install
