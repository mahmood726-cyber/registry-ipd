// engine.mjs — VERBATIM re-export of the repo's own engine (src/engine.js, UMD/CommonJS).
// No logic copied or edited; we import the published module unchanged and expose its API
// (reconstruct, reconstructEnsemble, classifyTier) plus the internals bag `_` (coxLogHR,
// kmFromIPD, medianFromKM, rmst, parametricArm, mulberry32, ...) for the harness.
import { createRequire } from 'node:module';
const require = createRequire(import.meta.url);
const RIPD = require('../src/engine.js');
export default RIPD;
export const { reconstruct, reconstructEnsemble, classifyTier } = RIPD;
