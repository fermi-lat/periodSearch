def generate(env, **kw):
    env.Tool('addLibrary', library = ['periodSearch'], package = 'periodSearch')
    env.Tool('pulsarDbLib')
    env.Tool('st_appLib')
    env.Tool('st_facilitiesLib')
    env.Tool('st_graphLib')
    env.Tool('st_streamLib')
    env.Tool('timeSystemLib')
    env.Tool('addLibrary', library = env['fftwLibs'])

def exists(env):
    return 1
