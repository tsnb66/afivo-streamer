OBJS := m_units_constants.o m_config.o m_lookup_table.o m_random.o		\
	m_photoi_mc.o m_streamer.o m_geometry.o m_transport_data.o m_field.o	\
	m_init_cond.o m_photoi_helmh.o m_photoi.o m_chemistry.o m_types.o m_gas.o

# Dependency information

m_chemistry.o: m_config.mod
m_chemistry.o: m_gas.mod
m_chemistry.o: m_lookup_table.mod
m_chemistry.o: m_types.mod
m_chemistry.o: m_units_constants.mod
m_field.o: m_chemistry.mod
m_field.o: m_streamer.mod
m_field.o: m_units_constants.mod
m_gas.o: m_config.mod
m_gas.o: m_types.mod
m_gas.o: m_units_constants.mod
m_init_cond.o: m_chemistry.mod
m_init_cond.o: m_geometry.mod
m_init_cond.o: m_streamer.mod
m_photoi_helmh.o: m_config.mod
m_photoi_helmh.o: m_gas.mod
m_photoi_helmh.o: m_streamer.mod
m_photoi_helmh.o: m_units_constants.mod
m_photoi_mc.o: m_config.mod
m_photoi_mc.o: m_gas.mod
m_photoi_mc.o: m_lookup_table.mod
m_photoi_mc.o: m_random.mod
m_photoi_mc.o: m_streamer.mod
m_photoi_mc.o: m_units_constants.mod
m_photoi.o: m_config.mod
m_photoi.o: m_gas.mod
m_photoi.o: m_photoi_helmh.mod
m_photoi.o: m_photoi_mc.mod
m_photoi.o: m_streamer.mod
m_photoi.o: m_units_constants.mod
m_streamer.o: m_chemistry.mod
m_streamer.o: m_config.mod
m_streamer.o: m_lookup_table.mod
m_streamer.o: m_random.mod
m_streamer.o: m_transport_data.mod
m_streamer.o: m_types.mod
m_streamer.o: m_units_constants.mod
streamer.o: m_chemistry.mod
streamer.o: m_field.mod
streamer.o: m_flux_schemes.mod
streamer.o: m_gas.mod
streamer.o: m_geometry.mod
streamer.o: m_init_cond.mod
streamer.o: m_photoi.mod
streamer.o: m_streamer.mod
streamer.o: m_units_constants.mod

# Hide some incorrect warnings
m_photoi_helmh.o: FFLAGS += -Wno-unused-function
m_photoi.o: FFLAGS += -Wno-unused-function

