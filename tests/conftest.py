from hypothesis import settings

settings.register_profile("ci", settings(max_examples=200))
settings.register_profile("dev", settings(max_examples=5))
