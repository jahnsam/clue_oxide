
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
class CluEOptions:
  def __init__(self):
    self._config = {}
    self._filters = {}
    self._spin_properties = {}
    self._structure_properties = {}
  #-----------------------------------------------------------------------------
  def set_config(self, key: str, value: str):
    self._config[key] = value
  #-----------------------------------------------------------------------------
  def add_filter(self, label: str):
    if label in self._filters:
      raise ValueError(f"filter[\"{label}\"] already exists")
    self._filters[label] = {}
  #-----------------------------------------------------------------------------
  def set_filter(self, label: str, key: str, value: str):
    if not label in self._filters:
      self.add_filter(label)
    self._filters[label][key] = value
  #-----------------------------------------------------------------------------
  def add_spin_properties(self, label: str, isotope):
    if label in self._spin_properties:
      if isotope in self._spin_properties[label]:
        raise ValueError(\
            f"spin_properties[\"{label}\"][\"{isotope}\"] already exists")
    self._spin_properties[label] = {isotope: {}}
  #-----------------------------------------------------------------------------
  def set_spin_property(self, label: str, isotope: str, key: str, value: str):
    if not label in self._spin_properties \
      or not isotope in self._spin_properties[label]:
      self.add_spin_properties(label,isotope)

    self._spin_properties[label][isotope][key] = value
  #-----------------------------------------------------------------------------
  def add_structure_properties(self, label: str):
    if label in self._structure_properties:
      raise ValueError(f"structure_properties[\"{label}\"] already exists")
    self._structure_properties[label] = {}
  #-----------------------------------------------------------------------------
  def set_structure_properties(self, label: str, key: str, value: str):
    if not label in self._structure_properties:
      self.add_structure_properties(label)
    self._structure_properties[label][key] = value
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

