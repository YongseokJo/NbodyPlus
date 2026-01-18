# Summary: 12-01 Restructure Main Loop for Async Pattern

## Commit
aad70b8

## Deliverables
- FEWBODY force loop restructured to use async MPI pattern
- Non-FEWBODY force loop also updated to async pattern
- Pre-posted receives via `postAllReceives()` before force loops
- Async sends via `runQueueAsync()` instead of blocking `runQueueAuto()`
- Non-blocking completion tests via `testQueueAsync()`
- Async callbacks via `callbackAsync()` inside inner loop

## Files Modified
- `src/irregular_routines.cpp` (lines 156-226, 250-259)

## Requirements Satisfied
- IRRG-01: Main loop restructured for async pattern
- IRRG-02: Multiple tasks sent before waiting for completions

## Technical Notes
- Callback moved inside inner loop to process completions as they arrive
- CM particle iteration happens during async window (overlap opportunity)
- Both FEWBODY and non-FEWBODY paths updated for consistency
