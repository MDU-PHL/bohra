"""
Automate deployment to PyPi
"""

import invoke


@invoke.task
def deploy(ctx):
    """
    Automate deployment
    rm -rf build/* dist/*
    bumpversion patch --verbose
    python3 setup.py sdist bdist_wheel
    twine upload dist/*
    git push --tags
    """
    ctx.run("rm -rf build/* dist/*")
    # ctx.run("bumpversion {bump} --verbose")
    ctx.run("python3 setup.py sdist bdist_wheel")
    ctx.run("python3 -m twine check dist/*")
    ctx.run("python3 -m twine upload dist/*")
    ctx.run("git push origin")
    ctx.run("git push kristy")

@invoke.task
def gitpush(ctx, message):
    """
    Automate deployment
    rm -rf build/* dist/*
    bumpversion patch --verbose
    python3 setup.py sdist bdist_wheel
    twine upload dist/*
    git push --tags
    """
    message = ' '.join(message.split('_'))
    ctx.run("git add -A")
    ctx.run(f"git commit -m '{message}'")
    ctx.run("git push origin")
    ctx.run("git push kristy")

