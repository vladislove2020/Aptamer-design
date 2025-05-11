"""
Это скрипт для сбора датасета из PDB для обучения генератора
"""

from typing import Dict, Any, Optional, Iterable
import asyncio
from tqdm.asyncio import tqdm

import requests
import json


# Пользуемся официальной апишкой PDB: https://data.rcsb.org/redoc/index.html
async def get_pdb_entry(pdb_id: str, assembly_id: int) -> Optional[Dict[str, Any]]:
  resp = requests.get(f"https://data.rcsb.org/rest/v1/core/assembly/{pdb_id}/{assembly_id}")
  if resp.status_code == 200:
      contents_str = resp.content.decode()
      try:
          contents_dict = json.loads(contents_str)
      except:
          contents_dict = None
      return contents_dict
  return None


async def get_polymers_info(pdb_id: str, polymer_id: int) -> Dict[str, Any]:
  resp = requests.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{polymer_id}")
  if resp.status_code == 200:
      contents_str = resp.content.decode()
      try:
          contents_dict = json.loads(contents_str)
      except:
          contents_dict = None
      return contents_dict
  return None


def get_all_pdb_ids():
    resp = requests.get("https://data.rcsb.org/rest/v1/holdings/current/entry_ids")
    contents_str = resp.content.decode()
    contents_dict = json.loads(contents_str)
    return contents_dict


async def get_protein_na_complexes(pdb_ids: Iterable[str]) -> Dict[str, Any]:
    result = {}

    async for pdb_id in tqdm(pdb_ids):
        entry_info = (await get_pdb_entry(pdb_id=pdb_id, assembly_id=1)) or {}

        polymer_entity_count = entry_info["rcsb_assembly_info"].get("polymer_entity_count", 0)
        polymer_entity_count_dna = entry_info["rcsb_assembly_info"].get("polymer_entity_count_dna", 0)
        polymer_entity_count_rna = entry_info["rcsb_assembly_info"].get("polymer_entity_count_rna", 0)
        polymer_entity_count_protein = entry_info["rcsb_assembly_info"].get("polymer_entity_count_protein", 0)
        if polymer_entity_count == 0:
            polymer_entity_count = polymer_entity_count_dna + polymer_entity_count_rna + polymer_entity_count_protein
        if polymer_entity_count == 0:
            continue

        if polymer_entity_count_dna == 0 and polymer_entity_count_rna == 0:
            continue

        new_line = {}

        new_line["pdb_id"] = pdb_id
        new_line["dna"] = (polymer_entity_count_dna > 0)
        new_line["rna"] = (polymer_entity_count_rna > 0)
        new_line["protein"] = (polymer_entity_count_protein > 0)

        chains_info = {}

        for entity_id in range(1, polymer_entity_count + 1):
            entity_info = (await get_polymers_info(pdb_id=pdb_id, polymer_id=entity_id)) or {}
            if len(entity_info.get("rcsb_polymer_entity_container_identifiers", {}).get("asym_ids", [])) != 1:
                continue
            chains_info[entity_info.get("rcsb_polymer_entity_container_identifiers", {}).get("asym_ids", [])[0]] = {
                "ids": entity_info.get("rcsb_polymer_entity_container_identifiers", {}).get("asym_ids", None),
                "auth_ids": entity_info.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids", None),
                "chem_comp_monomers": entity_info.get("rcsb_polymer_entity_container_identifiers", {}).get("chem_comp_monomers", None),
                "molecule_type": entity_info.get("entity_poly", {}).get("rcsb_entity_polymer_type", None),
                "sequence": entity_info.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", None),
                "entity_id": entity_id,
            }

        new_line["chains"] = chains_info

        result[pdb_id] = new_line

    return result


# Хоть бы не забанили за количество запросов (RPS можно подправить через asyncio.sleep мб)
async def main() -> None:
    total_pdb_ids = get_all_pdb_ids()

    result = await get_protein_na_complexes(pdb_ids=total_pdb_ids)

    with open("protein_na_complexes_info.json", "w") as js:
        json.dump(
            obj=result,
            fp=js,
            indent=4,
            ensure_ascii=False
        )
    
    return None


if __name__ == "__main__":
    asyncio.run(main())
