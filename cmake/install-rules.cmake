install(
    TARGETS NetworkV4_exe
    RUNTIME COMPONENT NetworkV4_Runtime
)

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif()
