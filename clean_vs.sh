rm -rf Debug Release
rm -rf glew SDL2
for SUBDIR in sdl2_opengl{22,33} ; do
    rm -rf ${SUBDIR}/Debug
    rm -rf ${SUBDIR}/Release
    rm -rf ${SUBDIR}/.vs
    rm -rf ${SUBDIR}/sdl_mandelbrot.vcxproj~
    rm -rf ${SUBDIR}/sdl_mandelbrot.vcxproj.user
    rm -rf ${SUBDIR}/.sdl_mandelbrot.vcxproj.un~
    rm -rf ${SUBDIR}/mandelbrot_gl.vcxproj~
    rm -rf ${SUBDIR}/mandelbrot_gl.vcxproj.user
    rm -rf ${SUBDIR}/.mandelbrot_gl.vcxproj.un~
done
