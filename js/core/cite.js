// Read an APA-style citation string into its parts, for labelling a study by its paper rather
// than by its serial number. Every field is null when the string does not have that shape;
// callers fall back to "Study N", never to a guess.
//   "Georgii, E., Kugler, K., ... (2019). The systems architecture ... . The Plant Cell, 31(2), ..."
/** The citation as published, without the curation notes some entries carry after it
 *  ("Set the title in italics...", "If the lab prefers personal attribution..."). Those notes
 *  were written for whoever formats the reference list and follow the first URL, so the
 *  reference ends there. Entries without a URL are returned whole. */
export function cleanCitation(citation) {
  if (!citation) return citation;
  const m = citation.match(/^[\s\S]*?https?:\/\/\S+?(?=[\s)]|$)/);
  return m ? m[0].replace(/[.,;]$/, '') : citation;
}

export function studyMeta(citation) {
  citation = cleanCitation(citation);
  const out = { title: null, firstAuthor: null, year: null, journal: null };
  if (!citation) return out;
  const y = citation.match(/\((\d{4})[a-z]?\)/);
  if (y) out.year = y[1];
  const surname = citation.split(',')[0].trim();
  if (surname && surname.length < 40) {
    const nAuthors = (citation.split('(')[0].match(/,\s*[A-Z]\./g) || []).length;
    out.firstAuthor = nAuthors > 1 ? `${surname} et al.` : surname;
  }
  const after = y ? citation.slice(citation.indexOf(y[0]) + y[0].length).replace(/^\.\s*/, '') : '';
  const m = after.match(/^(.+?[.?!])\s+([^.,]+?)[,.]/);
  if (m) {
    const t = m[1].replace(/\.$/, '');
    // An abbreviation inside a title ("No.") cuts it short, and a cut title has open brackets.
    if ((t.match(/\(/g) || []).length === (t.match(/\)/g) || []).length) out.title = t;
    out.journal = m[2].trim();
  }
  return out;
}
