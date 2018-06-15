unzip GFX_APIS/glew-2.1.0-win32.zip
mv glew-2.1.0 glew
unzip GFX_APIS/SDL2-devel-2.0.8-VC.zip
mv SDL2-2.0.8 SDL2
mkdir SDL2/include/SDL2
mv SDL2/include/* SDL2/include/SDL2 2>/dev/null
mkdir Debug Release
cp SDL2/lib/x86/SDL2.dll Debug
cp glew/bin/Release/Win32/glew32.dll Debug
cp SDL2/lib/x86/SDL2.dll Release
cp glew/bin/Release/Win32/glew32.dll Release
