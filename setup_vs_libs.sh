function EE() {
	echo $@
	$@
}
SLEEPAMT=0.717

EE unzip GFX_APIS/glew-2.1.0-win32.zip
sleep ${SLEEPAMT}
EE mv glew-2.1.0 glew
sleep ${SLEEPAMT}
EE unzip GFX_APIS/SDL2-devel-2.0.8-VC.zip
sleep ${SLEEPAMT}
EE mv SDL2-2.0.8 SDL2
sleep ${SLEEPAMT}
EE mkdir SDL2/include/SDL2
sleep ${SLEEPAMT}
EE mv SDL2/include/* SDL2/include/SDL2
sleep ${SLEEPAMT}
EE mkdir Debug Release
sleep ${SLEEPAMT}
EE cp SDL2/lib/x86/SDL2.dll Debug
sleep ${SLEEPAMT}
EE cp glew/bin/Release/Win32/glew32.dll Debug
sleep ${SLEEPAMT}
EE cp SDL2/lib/x86/SDL2.dll Release
sleep ${SLEEPAMT}
EE cp glew/bin/Release/Win32/glew32.dll Release
sleep ${SLEEPAMT}
EE cp -r Debug sdl2_opengl22
sleep ${SLEEPAMT}
EE cp -r Release sdl2_opengl22
sleep ${SLEEPAMT}
EE cp -r Debug sdl2_opengl33
sleep ${SLEEPAMT}
EE cp -r Release sdl2_opengl33

echo
echo "**** Development and Runtime libraries are now setup ****"
echo
