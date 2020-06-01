class NfUtils{
    def check_internal_overrides(String moduleName, Map params)
    {
        // get params set of keys
        Set paramsKeySet = params.keySet()

        // Interate through and set internals to the correct parameter at runtime
        paramsKeySet.each {
            if(it.startsWith("internal_")) {

                def searchString = moduleName + '_' + it.replace('internal_', '');

                if(paramsKeySet.contains(searchString)) {
                    params.replace(it, params.get(searchString))
                }
            }
        }
    }
}