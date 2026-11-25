// Boolean pipeline parameters.
//
// From the command line, `--skip_matam false` reaches the pipeline as the
// string "false" (a bare `--skip_matam` as "true"), and so does a quoted
// "false" in a params file. Groovy treats any non-empty string as true, so
// `if ( params.skip_matam )` would skip MATAM. Read boolean parameters with
// Flags.get(params, 'name') instead.
class Flags {
    static boolean get(Map params, String name) {
        def value = params[name]
        if ( value instanceof Boolean ) {
            return value
        }
        def text = value == null ? '' : value.toString().trim().toLowerCase()
        if ( text == 'true' ) {
            return true
        }
        if ( text == 'false' ) {
            return false
        }
        throw new IllegalArgumentException("--${name} must be true or false, not '${value}'")
    }
}
