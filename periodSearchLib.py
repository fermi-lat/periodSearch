#$Id$
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['periodSearch'])
    env.Tool('pulsarDbLib')
    env.Tool('st_appLib')
    env.Tool('st_facilitiesLib')
    env.Tool('st_graphLib')
    env.Tool('st_streamLib')
    env.Tool('timeSystemLib')
    env.Tool('addLibrary', library = env['fftwLibs'])

def exists(env):
    return 1
