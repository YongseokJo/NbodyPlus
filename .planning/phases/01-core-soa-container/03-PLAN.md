# Plan 03: Build Integration and Compilation Test

```yaml
wave: 3
depends_on: [01, 02]
files_modified:
  - src/Makefile
autonomous: true
```

## Objective

Integrate `particle_data.cpp` into the build system and verify compilation succeeds with no errors or warnings.

## Tasks

<task id="1">
Read current Makefile to understand build structure:
```bash
cat src/Makefile | head -100
```
Identify:
- Object file list variable
- Compilation rules
- Include paths
</task>

<task id="2">
Add `particle_data.o` to the object file list in Makefile.

Look for pattern like:
```makefile
OBJS = main.o default_global.o ...
```
Add `particle_data.o` to this list.
</task>

<task id="3">
Add compilation rule for particle_data.cpp (if not using pattern rule):
```makefile
particle_data.o: particle_data.cpp particle_data.h def.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
```
</task>

<task id="4">
Compile particle_data.cpp in isolation to check for errors:
```bash
cd src && g++ -std=c++11 -Wall -Wextra -c particle_data.cpp -o particle_data.o
```
Fix any compilation errors or warnings.
</task>

<task id="5">
Run full build to ensure no link errors:
```bash
cd src && make clean && make
```
(Note: At this phase, particle_data is not yet used by other files, so link should succeed)
</task>

<task id="6">
Verify no compiler warnings with strict flags:
```bash
g++ -std=c++11 -Wall -Wextra -Wpedantic -c src/particle_data.cpp -I src/ -o /dev/null
```
</task>

## Verification

- [ ] `particle_data.o` added to Makefile
- [ ] `particle_data.cpp` compiles without errors
- [ ] No compiler warnings with `-Wall -Wextra`
- [ ] Full build completes successfully

## must_haves

- [ ] ParticleData compiles standalone
- [ ] ParticleData integrated into Makefile
- [ ] No compilation warnings
