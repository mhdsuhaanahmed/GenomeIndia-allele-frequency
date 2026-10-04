def query_gnomad_variant(chrom, pos, ref, alt, max_retries=2):
    variant_id = f"{chrom.replace('chr','')}-{pos}-{ref}-{alt}"
    query = """
    query VariantQuery($variantId: String!, $datasetId: DatasetId!) {
      variant(variantId: $variantId, dataset: $datasetId) {
        genome { af }
        exome { af }
      }
    }
    """
    variables = {"variantId": variant_id, "datasetId": "gnomad_r4"}

    for attempt in range(max_retries):
        try:
            response = requests.post(
                "https://gnomad.broadinstitute.org/api",
                json={"query": query, "variables": variables},
                timeout=20
            )
            if response.status_code == 200:
                data = response.json()
                variant_data = data.get("data", {}).get("variant")
                if variant_data:
                    genome_af = variant_data.get("genome", {}).get("af") if variant_data.get("genome") else None
                    exome_af = variant_data.get("exome", {}).get("af") if variant_data.get("exome") else None
                    return genome_af, exome_af, "found_or_absent"
                else:
                    return None, None, "genuinely_absent"  # API responded but variant is null = truly not in gnomAD
        except requests.exceptions.RequestException as e:
            if attempt < max_retries - 1:
                time.sleep(1)
                continue
            return None, None, "query_failed"
    return None, None, "query_failed"