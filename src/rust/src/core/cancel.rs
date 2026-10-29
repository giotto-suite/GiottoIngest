//! Interrupting a long pass from R. The work runs on a worker thread while
//! R's main thread polls for a user interrupt (Ctrl-C, SIGINT); an interrupt
//! sets a flag the passes check between batches, tiles and grid cells, so
//! they stop early and the R side cleans up as after any error.

use std::ffi::c_void;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Duration;

static CANCEL: AtomicBool = AtomicBool::new(false);

/// `Err("interrupted")` once the user has interrupted the call.
pub fn check() -> Result<(), String> {
    if CANCEL.load(Ordering::Relaxed) { Err("interrupted".into()) } else { Ok(()) }
}

// R's API (Rinternals.h); not bound by extendr-ffi.
extern "C" {
    fn R_ToplevelExec(fun: extern "C" fn(*mut c_void), data: *mut c_void) -> i32;
    fn R_CheckUserInterrupt();
}

extern "C" fn check_interrupt(_: *mut c_void) {
    unsafe { R_CheckUserInterrupt() }
}

/// Whether R has a pending user interrupt. R_CheckUserInterrupt() jumps on
/// an interrupt; R_ToplevelExec() catches the jump and returns FALSE. Only
/// call this on R's main thread.
fn r_interrupted() -> bool {
    unsafe { R_ToplevelExec(check_interrupt, std::ptr::null_mut()) == 0 }
}

/// Runs `f` on a worker thread; this (R's main) thread polls for an
/// interrupt meanwhile. A panic in `f` becomes an error.
pub fn run<T: Send>(f: impl FnOnce() -> Result<T, String> + Send) -> Result<T, String> {
    CANCEL.store(false, Ordering::Relaxed);
    let out = std::thread::scope(|s| {
        let h = s.spawn(f);
        while !h.is_finished() {
            std::thread::sleep(Duration::from_millis(50));
            if !CANCEL.load(Ordering::Relaxed) && r_interrupted() {
                CANCEL.store(true, Ordering::Relaxed);
            }
        }
        h.join()
    });
    let interrupted = CANCEL.swap(false, Ordering::Relaxed);
    match out {
        Ok(r) if interrupted => r.and(Err("interrupted".into())),
        Ok(r) => r,
        Err(p) => Err(format!("internal error: {}", panic_message(&p))),
    }
}

fn panic_message(p: &Box<dyn std::any::Any + Send>) -> String {
    p.downcast_ref::<&str>().map(|s| s.to_string())
        .or_else(|| p.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| "panic".into())
}
