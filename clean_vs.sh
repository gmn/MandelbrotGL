
function EE() {
	echo $@
	$@
}

EE rm -rf Debug Release
EE rm -rf glew SDL2
for SUBDIR in sdl2_opengl{22,33} ; do
    EE rm -rf ${SUBDIR}/Debug
    EE rm -rf ${SUBDIR}/Release
    EE rm -rf ${SUBDIR}/.vs
    EE rm -rf ${SUBDIR}/sdl_mandelbrot.vcxproj~
    EE rm -rf ${SUBDIR}/sdl_mandelbrot.vcxproj.user
    EE rm -rf ${SUBDIR}/.sdl_mandelbrot.vcxproj.un~
    EE rm -rf ${SUBDIR}/mandelbrot_gl.vcxproj~
    EE rm -rf ${SUBDIR}/mandelbrot_gl.vcxproj.user
    EE rm -rf ${SUBDIR}/.mandelbrot_gl.vcxproj.un~
done
