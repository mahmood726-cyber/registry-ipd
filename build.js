/* Inline src/engine.js + examples/examples.js into a single fully-offline HTML file.
 * Produces dist/registry-ipd.html. Run: node build.js */
const fs = require('fs'), path = require('path');
const root = __dirname;
let html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
const engine = fs.readFileSync(path.join(root, 'src', 'engine.js'), 'utf8');
const examples = fs.readFileSync(path.join(root, 'examples', 'examples.js'), 'utf8');

// guard: never emit a literal </script> inside an inlined script block
function safe(js) { return js.replace(/<\/script>/g, '<\\/script>'); }

html = html.replace('<script src="src/engine.js"></script>', '<script>\n' + safe(engine) + '\n</script>');
html = html.replace('<script src="examples/examples.js"></script>', '<script>\n' + safe(examples) + '\n</script>');

if (html.includes('src="src/engine.js"') || html.includes('src="examples/')) {
  console.error('ERROR: external script refs remain — inlining failed'); process.exit(1);
}
fs.mkdirSync(path.join(root, 'dist'), { recursive: true });
const out = path.join(root, 'dist', 'registry-ipd.html');
fs.writeFileSync(out, html);
console.log('wrote', out, '(' + (html.length / 1024).toFixed(0) + ' KB, single-file offline)');
