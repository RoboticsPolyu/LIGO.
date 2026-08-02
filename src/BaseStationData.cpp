#include "BaseStationData.h"

#include <gnss_comm/gnss_constant.hpp>
#include <gnss_comm/gnss_utility.hpp>

#include <zip.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <utility>

namespace
{
namespace fs = std::filesystem;

std::string trim(std::string text)
{
  const auto first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return {};
  const auto last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

std::string lower(std::string text)
{
  std::transform(text.begin(), text.end(), text.begin(), [](unsigned char value) {
    return static_cast<char>(std::tolower(value));
  });
  return text;
}

uint32_t parseSatellite(const std::string &token)
{
  if (token.size() < 3) return 0;
  uint32_t system = 0;
  switch (token.front())
  {
    case 'G': system = SYS_GPS; break;
    case 'R': system = SYS_GLO; break;
    case 'E': system = SYS_GAL; break;
    case 'C': system = SYS_BDS; break;
    default: return 0;
  }
  return gnss_comm::sat_no(
      system, static_cast<uint32_t>(std::stoul(token.substr(1, 2))));
}

double carrierFrequency(char system, const std::string &observation_type,
                        uint32_t satellite,
                        const std::map<uint32_t, int> &glonass_channels)
{
  if (observation_type.size() != 3 || observation_type.front() != 'L')
    return 0.0;
  const char band = observation_type[1];
  if ((system == 'G' || system == 'E') && band == '1') return FREQ1;
  if (system == 'G' && band == '2') return FREQ2;
  if (system == 'E' && band == '7') return FREQ7;
  if ((system == 'G' || system == 'E') && band == '5') return FREQ5;
  if (system == 'E' && band == '6') return FREQ6;
  if (system == 'E' && band == '8') return FREQ8;
  // RINEX 3.02 uses band 1 for BeiDou B1I (1561.098 MHz). Some newer files
  // encode the same legacy signal as band 2, so accept both spellings.
  if (system == 'C' && (band == '1' || band == '2')) return FREQ1_BDS;
  if (system == 'C' && band == '7') return FREQ2_BDS;
  if (system == 'C' && band == '5') return FREQ5;
  if (system == 'C' && band == '6') return FREQ3_BDS;
  if (system == 'C' && band == '8') return FREQ8;
  if (system == 'R' && band == '1')
  {
    const auto channel = glonass_channels.find(satellite);
    if (channel == glonass_channels.end()) return 0.0;
    return FREQ1_GLO + channel->second * DFRQ1_GLO;
  }
  if (system == 'R' && band == '2')
  {
    const auto channel = glonass_channels.find(satellite);
    if (channel == glonass_channels.end()) return 0.0;
    return FREQ2_GLO + channel->second * DFRQ2_GLO;
  }
  if (system == 'R' && band == '3') return FREQ3_GLO;
  return 0.0;
}

bool isObservationMember(const std::string &name)
{
  const std::string value = lower(name);
  if (value.size() >= 4 && value.substr(value.size() - 4) == ".obs") return true;
  // Traditional RINEX names end in .YYo (for example hksc080q.25o).
  const auto dot = value.find_last_of('.');
  return dot != std::string::npos && value.size() >= dot + 4 &&
         value.back() == 'o';
}

struct RinexDestination
{
  Eigen::Vector3d position = Eigen::Vector3d::Zero();
  bool position_found = false;
  std::map<double, std::vector<BaseCarrierObservation>> epochs;
  std::set<std::string> time_systems;
};

bool parseRinexObservation(std::istream &input, const std::string &source,
                           double carrier_std_cycles,
                           RinexDestination &destination, std::string &error)
{
  std::map<char, std::vector<std::string>> observation_types;
  std::map<char, size_t> expected_type_counts;
  std::map<uint32_t, int> glonass_channels;
  char continued_system = 0;
  std::string time_system = "GPS";
  bool version_found = false;
  bool observation_file = false;
  bool header_complete = false;
  Eigen::Vector3d file_position = Eigen::Vector3d::Zero();
  bool file_position_found = false;
  size_t line_number = 0;
  std::string line;

  try
  {
    while (std::getline(input, line))
    {
      ++line_number;
      if (line.find("RINEX VERSION / TYPE") != std::string::npos)
      {
        const double version = std::stod(line.substr(0, 9));
        version_found = version >= 3.0 && version < 4.0;
        observation_file = line.find("OBSERVATION DATA") != std::string::npos;
      }
      else if (line.find("APPROX POSITION XYZ") != std::string::npos)
      {
        std::istringstream values(line.substr(0, std::min<size_t>(60, line.size())));
        values >> file_position.x() >> file_position.y() >> file_position.z();
        file_position_found = static_cast<bool>(values);
      }
      else if (line.find("SYS / # / OBS TYPES") != std::string::npos)
      {
        if (!line.empty() && line.front() != ' ')
        {
          continued_system = line.front();
          expected_type_counts[continued_system] =
              static_cast<size_t>(std::stoul(line.substr(3, 3)));
          observation_types[continued_system].clear();
        }
        if (continued_system == 0)
          throw std::runtime_error("observation-type continuation without a system");
        for (size_t offset = 7; offset + 3 <= std::min<size_t>(60, line.size());
             offset += 4)
        {
          const std::string type = trim(line.substr(offset, 3));
          if (!type.empty() && observation_types[continued_system].size() <
                                   expected_type_counts[continued_system])
            observation_types[continued_system].push_back(type);
        }
      }
      else if (line.find("GLONASS SLOT / FRQ #") != std::string::npos)
      {
        std::istringstream values(line.substr(0, std::min<size_t>(60, line.size())));
        std::string satellite_token;
        while (values >> satellite_token)
        {
          if (satellite_token.size() != 3 || satellite_token.front() != 'R')
            continue;
          int channel = 0;
          if (!(values >> channel)) break;
          const uint32_t satellite = parseSatellite(satellite_token);
          if (satellite != 0) glonass_channels[satellite] = channel;
        }
      }
      else if (line.find("TIME OF FIRST OBS") != std::string::npos)
      {
        std::istringstream values(line.substr(0, std::min<size_t>(60, line.size())));
        std::string token;
        while (values >> token) time_system = token;
      }
      else if (line.find("END OF HEADER") != std::string::npos)
      {
        header_complete = true;
        break;
      }
    }

    if (!version_found || !observation_file)
    {
      error = source + " is not a supported RINEX 3 observation file";
      return false;
    }
    if (!header_complete)
    {
      error = source + " has no END OF HEADER record";
      return false;
    }
    if (!file_position_found || file_position.norm() < 1.0)
    {
      error = source + " has no valid APPROX POSITION XYZ record";
      return false;
    }
    for (const auto &entry : expected_type_counts)
      if (observation_types[entry.first].size() != entry.second)
      {
        error = source + " has an incomplete SYS / # / OBS TYPES record";
        return false;
      }

    destination.time_systems.insert(time_system);

    if (destination.position_found &&
        (destination.position - file_position).norm() > 0.01)
    {
      error = source + " has a base position inconsistent with earlier files";
      return false;
    }
    destination.position = file_position;
    destination.position_found = true;

    std::string pending_line;
    bool has_pending_line = false;
    auto nextLine = [&](std::string &value) {
      if (has_pending_line)
      {
        value = std::move(pending_line);
        has_pending_line = false;
        return true;
      }
      if (!std::getline(input, value)) return false;
      ++line_number;
      return true;
    };

    while (nextLine(line))
    {
      if (trim(line).empty()) continue;
      if (line.front() != '>')
        throw std::runtime_error("expected an epoch record");

      // RINEX permits event records (for example flag 4 followed by header
      // updates) whose calendar fields are entirely blank.  Read the fixed
      // flag/count columns first so those records can be skipped without
      // requiring a timestamp.
      int event_flag = 0;
      int observation_count = 0;
      if (line.size() >= 35)
      {
        const std::string flag_text = trim(line.substr(31, 1));
        const std::string count_text = trim(line.substr(32, 3));
        if (!flag_text.empty()) event_flag = std::stoi(flag_text);
        if (!count_text.empty()) observation_count = std::stoi(count_text);
      }
      if (observation_count < 0)
        throw std::runtime_error("negative observation count");
      if (event_flag > 1)
      {
        for (int index = 0; index < observation_count; ++index)
          if (!nextLine(line)) throw std::runtime_error("incomplete event record");
        continue;
      }

      std::istringstream epoch_values(line.substr(1));
      double epoch[6] = {0};
      if (!(epoch_values >> epoch[0] >> epoch[1] >> epoch[2] >> epoch[3] >>
            epoch[4] >> epoch[5] >> event_flag >> observation_count))
        throw std::runtime_error("invalid epoch record");
      if (observation_count < 0)
        throw std::runtime_error("negative observation count");

      gnss_comm::gtime_t epoch_time = gnss_comm::epoch2time(epoch);
      if (time_system == "UTC") epoch_time = gnss_comm::utc2gpst(epoch_time);
      else if (time_system == "BDT") epoch_time = gnss_comm::time_add(epoch_time, 14.0);
      const double timestamp = gnss_comm::time2sec(epoch_time);
      std::map<std::pair<uint32_t, double>, BaseCarrierObservation> unique;

      for (int record_index = 0; record_index < observation_count; ++record_index)
      {
        if (!nextLine(line)) throw std::runtime_error("incomplete observation epoch");
        if (line.size() < 3) throw std::runtime_error("short observation record");
        const std::string satellite_token = line.substr(0, 3);
        const char system = satellite_token.front();
        const uint32_t satellite = parseSatellite(satellite_token);
        const auto types = observation_types.find(system);
        const size_t type_count =
            types == observation_types.end() ? 0 : types->second.size();
        std::string payload = line.size() > 3 ? line.substr(3) : std::string();

        // Most RINEX 3 writers use one long satellite line. Also accept strict
        // 80-column continuation lines, which begin with three blanks.
        while (payload.size() < type_count * 16)
        {
          std::string continuation;
          if (!nextLine(continuation)) break;
          if (continuation.size() >= 3 && continuation.substr(0, 3) == "   ")
            payload += continuation.substr(3);
          else
          {
            pending_line = std::move(continuation);
            has_pending_line = true;
            break;
          }
        }
        if (satellite == 0 || types == observation_types.end()) continue;

        for (size_t type_index = 0; type_index < type_count; ++type_index)
        {
          const std::string &type = types->second[type_index];
          const double frequency = carrierFrequency(
              system, type, satellite, glonass_channels);
          if (frequency <= 0.0) continue;
          const size_t offset = type_index * 16;
          if (offset >= payload.size()) continue;
          const std::string value_text =
              trim(payload.substr(offset, std::min<size_t>(14, payload.size() - offset)));
          if (value_text.empty()) continue;

          BaseCarrierObservation observation;
          observation.satellite = satellite;
          observation.frequency_hz = frequency;
          observation.carrier_cycles = std::stod(value_text);
          observation.carrier_std_cycles = carrier_std_cycles;
          if (offset + 14 < payload.size() &&
              std::isdigit(static_cast<unsigned char>(payload[offset + 14])))
            observation.loss_of_lock =
                static_cast<uint8_t>(payload[offset + 14] - '0');
          // Keep the first valid phase when a receiver publishes multiple
          // tracking codes on the same carrier frequency.
          unique.emplace(std::make_pair(satellite, frequency), observation);
        }
      }

      auto &destination_epoch = destination.epochs[timestamp];
      for (const auto &entry : unique) destination_epoch.push_back(entry.second);
    }
  }
  catch (const std::exception &exception)
  {
    error = source + ": parse error near line " + std::to_string(line_number) +
            ": " + exception.what();
    return false;
  }
  return true;
}

bool parseZipArchive(const std::string &path, double carrier_std_cycles,
                     RinexDestination &destination, std::string &error)
{
  int zip_error = 0;
  zip_t *archive = zip_open(path.c_str(), ZIP_RDONLY, &zip_error);
  if (!archive)
  {
    error = "cannot open ZIP archive " + path;
    return false;
  }
  std::unique_ptr<zip_t, decltype(&zip_close)> archive_guard(archive, zip_close);
  bool observation_found = false;
  const zip_int64_t entry_count = zip_get_num_entries(archive, 0);
  for (zip_uint64_t index = 0; index < static_cast<zip_uint64_t>(entry_count); ++index)
  {
    zip_stat_t status;
    zip_stat_init(&status);
    if (zip_stat_index(archive, index, 0, &status) != 0 || !status.name ||
        !isObservationMember(status.name))
      continue;
    zip_file_t *member = zip_fopen_index(archive, index, 0);
    if (!member)
    {
      error = "cannot read " + std::string(status.name) + " from " + path;
      return false;
    }
    std::unique_ptr<zip_file_t, decltype(&zip_fclose)> member_guard(member,
                                                                    zip_fclose);
    std::string content(static_cast<size_t>(status.size), '\0');
    zip_uint64_t offset = 0;
    while (offset < status.size)
    {
      const zip_int64_t count =
          zip_fread(member, &content[static_cast<size_t>(offset)], status.size - offset);
      if (count <= 0)
      {
        error = "failed while reading " + std::string(status.name);
        return false;
      }
      offset += static_cast<zip_uint64_t>(count);
    }
    std::istringstream input(content);
    if (!parseRinexObservation(input, path + ":" + status.name,
                               carrier_std_cycles, destination, error))
      return false;
    observation_found = true;
  }
  if (!observation_found)
  {
    error = path + " contains no RINEX observation member";
    return false;
  }
  return true;
}
}  // namespace

bool BaseStationData::load(const std::string &path, std::string &error,
                           double default_carrier_std_cycles)
{
  epochs_.clear();
  ecef_position_.setZero();
  error.clear();
  if (!(default_carrier_std_cycles > 0.0) ||
      !std::isfinite(default_carrier_std_cycles))
  {
    error = "default carrier standard deviation must be positive";
    return false;
  }

  RinexDestination destination;
  std::vector<fs::path> sources;
  const fs::path requested(path);
  if (!fs::exists(requested))
  {
    error = "cannot open " + path;
    return false;
  }
  if (fs::is_directory(requested))
  {
    for (const auto &entry : fs::directory_iterator(requested))
      if (entry.is_regular_file()) sources.push_back(entry.path());
    std::sort(sources.begin(), sources.end());
  }
  else
  {
    sources.push_back(requested);
  }

  bool parsed_any = false;
  for (const fs::path &source : sources)
  {
    const std::string extension = lower(source.extension().string());
    if (extension == ".zip")
    {
      if (!parseZipArchive(source.string(), default_carrier_std_cycles,
                           destination, error))
        return false;
      parsed_any = true;
    }
    else if (!fs::is_directory(requested) || isObservationMember(source.string()))
    {
      std::ifstream input(source);
      if (!input)
      {
        error = "cannot open " + source.string();
        return false;
      }
      if (!parseRinexObservation(input, source.string(),
                                 default_carrier_std_cycles, destination, error))
        return false;
      parsed_any = true;
    }
  }

  if (!parsed_any || !destination.position_found || destination.epochs.empty())
  {
    error = path + " contains no usable RINEX carrier observations";
    return false;
  }
  for (auto iterator = destination.epochs.begin();
       iterator != destination.epochs.end();)
  {
    if (iterator->second.empty()) iterator = destination.epochs.erase(iterator);
    else ++iterator;
  }
  if (destination.epochs.empty())
  {
    error = path + " contains no supported primary/secondary/L5 carrier observations";
    return false;
  }
  ecef_position_ = destination.position;
  epochs_ = std::move(destination.epochs);
  std::ostringstream time_systems;
  for (const std::string &system : destination.time_systems)
  {
    if (time_systems.tellp() > 0) time_systems << ',';
    time_systems << system;
  }
  time_system_ = time_systems.str().empty() ? "UNKNOWN" : time_systems.str();
  return true;
}

const std::vector<BaseCarrierObservation> *BaseStationData::epoch(
    double timestamp, double tolerance, double *matched_timestamp) const
{
  if (matched_timestamp) *matched_timestamp = 0.0;
  if (epochs_.empty()) return nullptr;
  auto after = epochs_.lower_bound(timestamp);
  auto best = after;
  if (after != epochs_.begin())
  {
    const auto before = std::prev(after);
    if (after == epochs_.end() ||
        std::abs(before->first - timestamp) <= std::abs(after->first - timestamp))
      best = before;
  }
  if (best == epochs_.end()) return nullptr;
  if (matched_timestamp) *matched_timestamp = best->first;
  if (std::abs(best->first - timestamp) > tolerance)
    return nullptr;
  return &best->second;
}

double BaseStationData::firstTimestamp() const
{
  return epochs_.empty() ? 0.0 : epochs_.begin()->first;
}

double BaseStationData::lastTimestamp() const
{
  return epochs_.empty() ? 0.0 : epochs_.rbegin()->first;
}
